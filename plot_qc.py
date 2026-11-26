"""Draw teaching figures from a completed pilot workspace or archive replay."""

import argparse
import csv
from pathlib import Path

import nibabel as nib
import numpy as np

from reproduce import on_grid, plt, values


def slice_index(mask, fraction):
    indices = np.flatnonzero(mask.any(axis=(0, 1)))
    if not len(indices):
        raise ValueError("The display mask is empty")
    return int(indices[round(fraction * (len(indices) - 1))])


def brain_slice(ax, image, mask, index, spacing):
    # Restrict the anatomical display to the brain mask.
    view = np.where(mask[:, :, index], image[:, :, index], np.nan)
    upper = np.percentile(image[mask], 98)
    ax.imshow(view.T, origin="lower", cmap="gray", vmin=0, vmax=upper,
              interpolation="nearest", aspect=spacing[1] / spacing[0])
    xx, yy = np.where(mask[:, :, index])
    ax.set(xlim=(xx.min() - 6, xx.max() + 6), ylim=(yy.min() - 6, yy.max() + 6),
           xticks=[], yticks=[])
    ax.set_facecolor("#101b26")


def contour(ax, mask, index, colour):
    plane = mask[:, :, index].T
    if plane.any() and not plane.all():
        ax.contour(plane, levels=[0.5], colors=[colour], linewidths=0.85)


def masks(work, out):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.8), layout="constrained")
    for posture, row in zip(("supine", "upright"), axes):
        folder = work / posture
        t1 = nib.load(folder / "t1_brain.nii.gz")
        t1_mask = on_grid(folder / "t1_brain_mask.nii.gz", t1) > 0
        t1_z = slice_index(t1_mask, 0.52)
        brain_slice(row[0], values(t1), t1_mask, t1_z, t1.header.get_zooms())
        contour(row[0], t1_mask, t1_z, "#65e1de")
        row[0].set_title(f"{posture.title()} · T1 BET · z {t1_z + 1}")
        mean = nib.load(folder / "mean.nii.gz")
        epi_mask = on_grid(folder / "epi_brain_mask.nii.gz", mean) > 0
        signal = values(mean)
        for ax, fraction in zip(row[1:], (0.26, 0.62)):
            z = slice_index(epi_mask, fraction)
            brain_slice(ax, signal, epi_mask, z, mean.header.get_zooms())
            contour(ax, epi_mask, z, "#65e1de")
            ax.set_title(f"BOLD BET · z {z + 1}")
    fig.suptitle("Volunteer 1 | T1 and mean BOLD brain masks", fontsize=16)
    fig.get_layout_engine().set(h_pad=0.06)
    fig.savefig(out / "09-brain-masks.png", dpi=160)
    plt.close(fig)


def atlas_overlay(work, out):
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.3), layout="constrained")
    for posture, row in zip(("supine", "upright"), axes):
        folder = work / posture
        ref = nib.load(folder / "mean.nii.gz")
        mean = values(ref)
        brain = on_grid(folder / "trim3.ica" / "mask.nii.gz", ref) > 0
        atlas = on_grid(folder / "atlas_native.nii.gz", ref).astype(int)
        usable = brain & (mean >= 0.25 * np.median(mean[brain]))
        index = slice_index(brain, 0.36)
        for ax in row:
            brain_slice(ax, mean, brain, index, ref.header.get_zooms())
        outline = brain & (atlas > 0)
        labels = np.ma.masked_where(~outline[:, :, index].T, atlas[:, :, index].T % 20)
        row[1].imshow(labels, origin="lower", cmap="tab20", vmin=0, vmax=19,
                      alpha=0.55, interpolation="nearest", aspect=ref.header.get_zooms()[1] / ref.header.get_zooms()[0])
        contour(row[2], brain, index, "#f3bc75")
        contour(row[2], usable, index, "#65e1de")
        row[0].set_title(f"{posture.title()} · mean BOLD · z {index + 1}")
        row[1].set_title("Atlas labels on native grid")
        row[2].set_title("ICA mask and usable signal")
    fig.suptitle("Volunteer 1 | Native-space atlas and signal coverage\n"
                 "Amber = ICA mask  ·  Cyan = signal at least 25% of brain median", fontsize=13)
    fig.savefig(out / "10-atlas-and-signal.png", dpi=160)
    plt.close(fig)


def parcel_coverage(work, out):
    with (work / "qc" / "parcels.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    supine = np.array([float(row["supine_fraction"]) for row in rows])
    upright = np.array([float(row["upright_fraction"]) for row in rows])
    kept = np.array([row["in_pair"] == "True" for row in rows])
    fig, ax = plt.subplots(figsize=(6.3, 5.8), layout="constrained")
    ax.axvspan(0.7, 1.04, ymin=(0.7 / 1.04), color="#e4f2ef", zorder=0)
    ax.axvline(0.7, color="0.5", linestyle=":")
    ax.axhline(0.7, color="0.5", linestyle=":")
    ax.scatter(supine[~kept], upright[~kept], s=27, color="#a94b2e", label="Excluded")
    ax.scatter(supine[kept], upright[kept], s=26, color="#216583", label="Retained in both")
    for i, row in enumerate(rows):
        if row["region"] == "Frontal Pole":
            ax.annotate("Frontal pole", (supine[i], upright[i]), xytext=(9, -24),
                        textcoords="offset points", arrowprops={"arrowstyle": "-", "color": "0.4"})
    ax.set(xlim=(0, 1.04), ylim=(0, 1.04), aspect="equal",
           xlabel="Supine usable atlas fraction", ylabel="Upright usable atlas fraction",
           title=f"Volunteer 1 | {kept.sum()} of {len(rows)} parcels shared")
    ax.legend(loc="lower right", title="Also require 50 native voxels", fontsize=8)
    fig.savefig(out / "11-parcel-screen.png", dpi=160)
    plt.close(fig)


def split_matrices(work, out):
    fig, axes = plt.subplots(2, 4, figsize=(11.5, 6.2), layout="constrained")
    order = (("before_wm", "first"), ("before_wm", "second"),
             ("after_wm", "first"), ("after_wm", "second"))
    for posture, row in zip(("supine", "upright"), axes):
        for ax, (arm, half) in zip(row, order):
            matrix = np.loadtxt(work / "qc" / f"{posture}_{arm}_{half}.csv",
                                delimiter=",", skiprows=1, usecols=None, dtype=str, quotechar='"')
            values_r = matrix[:, 1:].astype(float)
            np.fill_diagonal(values_r, np.nan)
            panel = ax.imshow(values_r, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
            ax.set(title=f"{'Before' if arm == 'before_wm' else 'After'} WM · {half} half",
                   xticks=[], yticks=[])
        row[0].set_ylabel(posture.title())
    fig.colorbar(panel, ax=axes, shrink=0.74, label="Pearson r")
    fig.suptitle("Volunteer 1 | The same shared parcels in every split-half matrix\n"
                 "All panels use ICA-cleaned BOLD. The diagonal is hidden.", fontsize=14)
    fig.savefig(out / "12-split-half-matrices.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    out = args.work / "qc" / "teaching_figures"
    out.mkdir(parents=True, exist_ok=True)
    for draw in (masks, atlas_overlay, parcel_coverage, split_matrices):
        draw(args.work, out)


if __name__ == "__main__":
    main()
