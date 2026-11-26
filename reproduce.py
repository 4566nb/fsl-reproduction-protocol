"""Reproduce the FSL pilot QC and split-half connectivity from local derivatives."""

import argparse
import csv
from pathlib import Path
from xml.etree import ElementTree

import matplotlib
import nibabel as nib
import numpy as np
from scipy.ndimage import (
    affine_transform,
    binary_closing,
    binary_erosion,
    binary_fill_holes,
    distance_transform_edt,
    gaussian_filter,
    label,
)
from scipy.stats import spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TR = 3.8
FRAMES = np.arange(4, 251)
HALVES = (slice(0, 123), slice(124, 247))


def values(image):
    return image.get_fdata(dtype=np.float32, caching="unchanged")


def on_grid(path, reference):
    image = nib.load(path)
    if image.shape[:3] != reference.shape[:3] or not np.allclose(
        image.affine, reference.affine
    ):
        raise ValueError(f"Grid mismatch: {path}")
    return values(image)


def write_table(path, header, rows, delimiter=","):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter=delimiter)
        writer.writerow(header)
        writer.writerows(rows)


def fsl_coordinates(image):
    """Map NIfTI voxel indices to the scaled voxel coordinates used by FLIRT."""
    spacing = image.header.get_zooms()[:3]
    transform = np.diag([*spacing, 1.0])
    if np.linalg.det(image.affine[:3, :3]) > 0:
        transform[0, 0] *= -1
        transform[0, 3] = (image.shape[0] - 1) * spacing[0]
    return transform


def bridge(args):
    native = nib.load(args.native)
    isotropic = nib.load(args.isotropic)
    native_to_iso = (
        fsl_coordinates(isotropic)
        @ np.linalg.inv(isotropic.affine)
        @ native.affine
        @ np.linalg.inv(fsl_coordinates(native))
    )
    np.savetxt(args.output, np.loadtxt(args.iso_to_mni) @ native_to_iso)


def motion(args):
    args.out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, layout="constrained")
    rows = []
    for posture, path, ax in zip(("supine", "upright"), (args.supine, args.upright), axes):
        parameters = np.loadtxt(path)
        if parameters.shape != (250, 6) or not np.isfinite(parameters).all():
            raise ValueError(f"Expected 250 finite MCFLIRT rows with six columns: {path}")
        delta = np.abs(np.diff(parameters, axis=0))
        fd = delta[:, 3:].sum(axis=1) + 50 * delta[:, :3].sum(axis=1)
        time = np.arange(1, 250) * TR / 60
        write_table(
            args.out / f"{posture}_fd.tsv",
            ("original_frame", "minutes_from_start", "power_fd_mm"),
            zip(range(1, 251), np.r_[0, time], np.r_[np.nan, fd]),
            delimiter="\t",
        )
        rows.append((posture, len(fd), fd.mean(), np.count_nonzero(fd > 0.2)))
        ax.plot(time, fd)
        ax.axhline(0.2, color="0.5", linestyle=":")
        ax.set(title=posture, ylabel="Power FD (mm)", ylim=(0, None))
    axes[-1].set_xlabel("Minutes from scan start")
    fig.savefig(args.out / "motion.png", dpi=160)
    plt.close(fig)
    write_table(args.out / "motion.csv", ("posture", "transitions", "mean_fd_mm", "above_0p2"), rows)


def ica_summary(path, out):
    image = nib.load(path / "melodic_IC.nii.gz")
    maps = values(image)
    mix = np.loadtxt(path / "melodic_mix", ndmin=2)
    if mix.shape[1] != maps.shape[3]:
        raise ValueError(f"Mixing matrix does not match component maps: {path}")
    spatial_energy = (maps**2).sum(axis=(0, 1), dtype=np.float64)
    spatial_fraction = spatial_energy.max(axis=0) / spatial_energy.sum(axis=0)
    time_energy = mix**2
    starts = time_energy[:3].sum(axis=0) / time_energy.sum(axis=0)
    write_table(
        out,
        ("component", "native_slice", "max_slice_energy_fraction", "first_three_energy_fraction"),
        zip(range(1, mix.shape[1] + 1), spatial_energy.argmax(axis=0) + 1, spatial_fraction, starts),
    )
    print(f"{path.name}: {mix.shape[0]} volumes, {mix.shape[1]} ICs, "
          f"{np.count_nonzero(spatial_fraction > 0.5)} single-slice flags, "
          f"max early energy {starts.max():.1%}")
    return image, maps


def ica(args):
    args.out.mkdir(parents=True, exist_ok=True)
    initial, before = ica_summary(args.initial, args.out / "ica_initial.csv")
    trimmed, after = ica_summary(args.trimmed, args.out / "ica_trimmed.csv")
    if initial.shape[:3] != trimmed.shape[:3] or not np.allclose(initial.affine, trimmed.affine):
        raise ValueError("Initial and trimmed ICA maps must be on the same grid")
    mask = (on_grid(args.initial / "mask.nii.gz", initial) > 0) & (
        on_grid(args.trimmed / "mask.nii.gz", initial) > 0
    )
    correlations = np.abs(np.corrcoef(before[mask].T, after[mask].T)[:before.shape[3], before.shape[3]:])
    write_table(
        args.out / "ica_spatial_matches.csv",
        ("initial_ic", "closest_trimmed_ic", "absolute_spatial_r"),
        zip(range(1, before.shape[3] + 1), correlations.argmax(axis=1) + 1, correlations.max(axis=1)),
    )


def slice_qc(args):
    args.out.mkdir(parents=True, exist_ok=True)
    reference = nib.load(args.mask)
    brain = values(reference) > 0
    sds = []
    for path in (args.before, args.after):
        data = on_grid(path, reference)
        if data.shape != (*brain.shape, len(FRAMES)):
            raise ValueError(f"Expected 247 volumes: {path}")
        sd = np.full(brain.shape[2], np.nan)
        for z in range(brain.shape[2]):
            core = binary_erosion(brain[:, :, z], iterations=3)
            voxels = data[:, :, z, :][core]
            if len(voxels) < 100:
                continue
            baseline = voxels.mean(axis=1)
            valid = baseline > 10
            if valid.any():
                percent = 100 * (voxels[valid] / baseline[valid, None] - 1)
                sd[z] = np.median(percent, axis=0).std()
        sds.append(sd)
    reduction = np.full_like(sds[0], np.nan)
    np.divide(sds[1], sds[0], out=reduction, where=sds[0] > 0)
    reduction = 100 * (1 - reduction)
    slices = np.arange(1, brain.shape[2] + 1)
    write_table(args.out / "slice_sd.csv",
                ("native_slice", "before_sd_pct", "after_sd_pct", "relative_reduction_pct"),
                zip(slices, *sds, reduction))
    fig, ax = plt.subplots(layout="constrained")
    ax.plot(slices, sds[0], "o-", label="Before IC regression", markersize=3)
    ax.plot(slices, sds[1], "o-", label="After IC regression", markersize=3)
    ax.set(xlabel="Native slice (1-based)", ylabel="Temporal SD (%)")
    ax.legend()
    fig.savefig(args.out / "slice_sd.png", dpi=160)
    plt.close(fig)


def atlas_names(path):
    return {
        int(node.attrib["index"]) + 1: node.text.strip()
        for node in ElementTree.parse(path).findall(".//label")
    }


def deep_white_matter(wm, usable, spacing):
    return usable & (wm >= 0.95) & (distance_transform_edt(wm >= 0.80, sampling=spacing) >= 3)


def equal_slice_trace(data, mask, slices):
    traces = [
        data[:, :, z, :][mask[:, :, z]].mean(axis=0, dtype=np.float64)
        for z in slices
    ]
    return np.mean(traces, axis=0)


def raw_tissue(raw, brain, native_atlas, wm, accepted, spacing, out):
    """Form the supine raw tissue masks and use identical slice support in each."""
    mean = raw.mean(axis=3)
    median_brain = np.median(mean[brain])
    cutoff = 0.25 * median_brain
    usable = brain & (mean >= cutoff)
    cortex = usable & (wm < 0.25) & np.isin(native_atlas, accepted)
    deep = deep_white_matter(wm, usable, spacing)
    background = np.zeros_like(brain)
    x, y = np.indices(brain.shape[:2])
    for z in range(5, min(23, brain.shape[2])):
        in_brain = brain[:, :, z]
        if not in_brain.any():
            continue
        plane = mean[:, :, z]
        head = binary_fill_holes(binary_closing(gaussian_filter(plane, 1) > cutoff, iterations=4))
        clusters, count = label(head)
        overlap = np.bincount(clusters[in_brain], minlength=count + 1)
        overlap[0] = 0
        if overlap.max() == 0:
            continue
        head = binary_fill_holes(clusters == overlap.argmax()) | in_brain
        distance = distance_transform_edt(~head, sampling=spacing[:2])
        background[:, :, z] = (
            (distance >= 12) & (x >= 4) & (x < brain.shape[0] - 4)
            & (y >= 4) & (y < brain.shape[1] - 4) & (plane >= 1) & (plane <= 0.18 * median_brain)
        )
    slices = np.flatnonzero(
        (cortex.sum(axis=(0, 1)) >= 100)
        & (deep.sum(axis=(0, 1)) >= 50)
        & (background.sum(axis=(0, 1)) >= 100)
    )
    if len(slices) < 8:
        raise ValueError("Supine raw tissue QC has fewer than eight shared slices")
    deep[:, :, ~np.isin(np.arange(brain.shape[2]), slices)] = False
    time = (FRAMES - 1) * TR / 60
    fig, axes = plt.subplots(3, 1, sharex=True, layout="constrained", figsize=(8, 5))
    rows = []
    for mask, name, ax in zip((cortex, deep, background), ("Cortex", "Deep WM", "Outside head"), axes):
        trace = equal_slice_trace(raw, mask, slices)
        percent = 100 * (trace / trace.mean() - 1)
        rows.append((name, len(slices), *(percent[half].std() for half in HALVES)))
        ax.plot(time, percent, linewidth=0.8)
        ax.axvspan(8.04, 9.94, color="#78aec3", alpha=0.2)
        ax.set(title=name, ylabel="Change (%)")
    axes[-1].set_xlabel("Minutes from scan start")
    fig.savefig(out / "supine_raw_tissue.png", dpi=160)
    plt.close(fig)
    write_table(out / "supine_raw_tissue.csv", ("tissue", "matched_slices", "first_sd_pct", "second_sd_pct"), rows)

    # Show the sampling regions while limiting anatomical intensities to the brain.
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.5), layout="constrained")
    colours = matplotlib.colors.ListedColormap(["#38c9de", "#f5a24b", "#ab76cf"])
    tissue = cortex.astype(np.uint8) + 2 * deep.astype(np.uint8) + 3 * background.astype(np.uint8)
    upper = np.percentile(mean[brain], 98)
    for ax, target in zip(axes, (8, 14, 20)):
        z = slices[np.argmin(abs(slices - target))]
        grey = np.where(brain[:, :, z], mean[:, :, z], np.nan)
        ax.imshow(grey.T, origin="lower", cmap="gray", vmin=0, vmax=upper)
        ax.set_facecolor("#17212b")
        overlay = np.ma.masked_equal(tissue[:, :, z].T, 0)
        ax.imshow(overlay, origin="lower", cmap=colours, vmin=0.5, vmax=3.5,
                  alpha=0.7, interpolation="nearest")
        ax.set(title=f"Native slice {z + 1}", xticks=[], yticks=[])
    fig.suptitle("Supine raw BOLD | Cortex (cyan), deep WM (orange), outside-head background (purple)")
    fig.savefig(out / "raw_tissue_masks.png", dpi=160)
    plt.close(fig)
    return deep


def atlas_lookup(ref, standard, mni_to_bold):
    return np.linalg.inv(fsl_coordinates(ref)) @ mni_to_bold @ fsl_coordinates(standard)


def nearest(data, lookup, shape):
    return affine_transform(data, lookup, output_shape=shape, mode="grid-constant", order=0)


def extract_run(folder, posture, standard, standard_labels, names, out):
    ref = nib.load(folder / "mean.nii.gz")
    mean = values(ref)
    brain = on_grid(folder / "trim3.ica" / "mask.nii.gz", ref) > 0
    atlas = on_grid(folder / "atlas_native.nii.gz", ref)
    wm = on_grid(folder / "wm_native.nii.gz", ref)
    if not brain.any() or not set(np.unique(atlas)) <= ({0} | names.keys()):
        raise ValueError(f"Check brain mask and atlas label IDs: {folder}")
    usable = brain & (mean >= 0.25 * np.median(mean[brain]))
    lookup = atlas_lookup(ref, standard, np.loadtxt(folder / "mni_to_bold.mat"))
    in_mni = nearest(usable, lookup, standard.shape)
    coverage = {i: in_mni[standard_labels == i].mean() for i in names}
    counts = {i: np.count_nonzero(usable & (atlas == i)) for i in names}
    accepted = [i for i in names if coverage[i] >= 0.70 and counts[i] >= 50]
    if len(accepted) < 2:
        raise ValueError(f"Fewer than two usable parcels: {folder}")
    spacing = ref.header.get_zooms()[:3]
    if posture == "supine":
        raw = on_grid(folder / "bold_raw.nii.gz", ref)
        if raw.shape != (*mean.shape, 250):
            raise ValueError("Expected 250 raw supine volumes")
        wm_mask = raw_tissue(raw[..., 3:], brain, atlas, wm, accepted, spacing, out)
        del raw
    else:
        wm_mask = deep_white_matter(wm, usable, spacing)
        wm_mask[:, :, :5] = False
        wm_mask[:, :, 22:] = False
    wm_mask &= ~np.isin(atlas, accepted)
    slices = np.flatnonzero(wm_mask.sum(axis=(0, 1)) >= 50)
    if len(slices) < 8:
        raise ValueError(f"Fewer than eight eligible white-matter slices: {folder}")
    data = on_grid(folder / "clean.nii.gz", ref)
    if data.shape != (*mean.shape, len(FRAMES)):
        raise ValueError(f"Expected 247 cleaned volumes: {folder}")
    wm_trace = equal_slice_trace(data, wm_mask, slices)
    nodes = {
        i: data[usable & (atlas == i)].mean(axis=0, dtype=np.float64)
        for i in accepted
    }
    return coverage, counts, nodes, wm_trace


def connectivity(args):
    args.out.mkdir(parents=True, exist_ok=True)
    standard = nib.load(args.atlas)
    standard_labels = values(standard)
    names = atlas_names(args.labels)
    if not set(names) <= set(np.unique(standard_labels)):
        raise ValueError("The atlas XML and standard maximum-probability label image must match")
    supine = extract_run(args.supine, "supine", standard, standard_labels, names, args.out)
    upright = extract_run(args.upright, "upright", standard, standard_labels, names, args.out)
    common = [i for i in names if i in supine[2] and i in upright[2]]
    if len(common) < 3:
        raise ValueError("At least three common parcels are required for edge-pattern correlation")
    parcel_names = [names[i] for i in common]
    write_table(
        args.out / "parcels.csv",
        ("atlas_label", "region", "supine_fraction", "supine_voxels",
         "upright_fraction", "upright_voxels", "in_pair"),
        ((i, name, supine[0][i], supine[1][i], upright[0][i], upright[1][i], i in common)
         for i, name in names.items()),
    )
    edges = np.triu_indices(len(common), k=1)
    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(8, 7), layout="constrained")
    for (posture, (_, _, traces, wm)), panels in zip(
        (("supine", supine), ("upright", upright)), axes
    ):
        nodes = np.column_stack([traces[i] for i in common])
        centered_wm = wm - wm.mean()
        if not np.isfinite(wm).all() or centered_wm @ centered_wm == 0:
            raise ValueError(f"Invalid deep WM time series: {posture}")
        slopes = centered_wm @ (nodes - nodes.mean(axis=0)) / (centered_wm @ centered_wm)
        adjusted = nodes - np.outer(centered_wm, slopes)
        for (condition, series), ax in zip((("before_wm", nodes), ("after_wm", adjusted)), panels):
            write_table(
                args.out / f"{posture}_{condition}_nodes.tsv",
                ("original_frame", "deep_wm_mean", *parcel_names),
                zip(FRAMES, wm, *series.T), delimiter="\t",
            )
            edge_vectors = []
            for half_name, half in zip(("first", "second"), HALVES):
                if not np.isfinite(series[half]).all() or np.any(series[half].std(axis=0) == 0):
                    raise ValueError(f"Invalid parcel series: {posture}, {condition}, {half_name}")
                matrix = np.corrcoef(series[half], rowvar=False)
                edge_vectors.append(matrix[edges])
                write_table(
                    args.out / f"{posture}_{condition}_{half_name}.csv",
                    ("parcel", *parcel_names),
                    ((name, *row) for name, row in zip(parcel_names, matrix)),
                )
            a, b = edge_vectors
            r = np.corrcoef(a, b)[0, 1]
            rows.append((posture, condition, len(common), len(a), np.median(a), np.median(b),
                         np.abs(a - b).mean(), r, spearmanr(a, b).statistic))
            ax.scatter(a, b, s=9, alpha=0.4)
            ax.plot([-1, 1], [-1, 1], ":", color="0.5")
            ax.set(xlim=(-1, 1), ylim=(-1, 1), aspect="equal", xlabel="First-half edge r",
                   ylabel="Second-half edge r", title=f"{posture}, {condition.replace('_', ' ')} · r={r:.3f}")
    write_table(
        args.out / "connectivity.csv",
        ("posture", "condition", "parcels", "edges", "first_median_r", "second_median_r",
         "edge_mae", "pattern_pearson_r", "pattern_spearman_r"),
        rows,
    )
    fig.savefig(args.out / "split_half.png", dpi=160)
    plt.close(fig)
    print(f"Accepted parcels: supine {len(supine[2])}, upright {len(upright[2])}, shared {len(common)}")


def coverage(args):
    """Screen Volunteer 2's raw BOLD coverage against Volunteer 1's paired nodes."""
    atlas = nib.load(args.atlas)
    standard = values(atlas).astype(int)
    with args.pilot_parcels.open(encoding="utf-8-sig", newline="") as stream:
        parcels = list(csv.DictReader(stream))
    label_ids = np.array([int(parcel["atlas_label"]) for parcel in parcels])
    header = ["atlas_label", "region", "v1_supine_fraction", "v1_upright_fraction", "v1_in_pair"]
    denominator = np.bincount(standard.ravel())
    measurements = []
    for posture, folder in (("supine", args.supine), ("upright", args.upright)):
        ref = nib.load(folder / "mean.nii.gz")
        mean = values(ref)
        brain = on_grid(folder / "epi_brain_mask.nii.gz", ref) > 0
        if not np.isfinite(mean).all() or not brain.any():
            raise ValueError(f"Invalid mean BOLD or mask: {folder}")
        usable = brain & (mean > 0.25 * np.median(mean[brain]))
        lookup = atlas_lookup(ref, atlas, np.loadtxt(folder / "mni_to_bold.mat"))
        usable_mni = nearest(usable, lookup, atlas.shape)
        fov_mni = nearest(np.ones(ref.shape, dtype=bool), lookup, atlas.shape)
        native_labels = nearest(standard, np.linalg.inv(lookup), ref.shape)
        usable_count = np.bincount(native_labels[usable], minlength=len(denominator))
        usable_mni_count = np.bincount(standard[usable_mni], minlength=len(denominator))
        fov_count = np.bincount(standard[fov_mni], minlength=len(denominator))
        fractions = usable_mni_count[label_ids] / denominator[label_ids]
        counts = usable_count[label_ids]
        passes = (fractions >= 0.70) & (counts >= 50)
        measurements.append((fov_count[label_ids] / denominator[label_ids], fractions, counts, passes))
        header.extend(f"v2_{posture}_{key}" for key in
                      ("fov_fraction", "usable_fraction", "usable_voxels", "passes"))

    supine, upright = measurements
    v1_pair = np.array([p["in_pair"] == "True" for p in parcels])
    v2_pair = supine[-1] & upright[-1]
    all_four = v1_pair & v2_pair
    header.append("candidate_all_four")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_table(
        args.out, header,
        ((label_id, p["region"], p["supine_fraction"], p["upright_fraction"], in_v1,
          *s, *u, in_all)
         for label_id, p, in_v1, s, u, in_all in zip(
             label_ids, parcels, v1_pair, zip(*supine), zip(*upright), all_four
         )),
    )
    print(f"V2 coverage passes: supine {supine[-1].sum()}, upright {upright[-1].sum()}, "
          f"both {v2_pair.sum()}. Common to all four scans: {all_four.sum()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    for command, function, paths in (
        ("bridge", bridge, ("native", "isotropic", "iso-to-mni", "output")),
        ("motion", motion, ("supine", "upright", "out")),
        ("ica", ica, ("initial", "trimmed", "out")),
        ("slice-qc", slice_qc, ("before", "after", "mask", "out")),
        ("connectivity", connectivity, ("supine", "upright", "atlas", "labels", "out")),
        ("coverage", coverage, ("supine", "upright", "atlas", "pilot-parcels", "out")),
    ):
        options = sub.add_parser(command)
        for name in paths:
            options.add_argument(f"--{name}", type=Path, required=True)
        options.set_defaults(run=function)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
