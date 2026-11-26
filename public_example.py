"""Generate the public protocol figures from artificial inputs only."""

from argparse import ArgumentParser
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import map_coordinates

from reproduce import HALVES, TR, motion, plt


BLUE = "#246982"
ORANGE = "#be6333"
CYAN = "#3dd5d3"
RNG = np.random.default_rng(618)


def save(fig, out, name, title):
    fig.suptitle(f"{title}\nSynthetic teaching example", fontsize=14)
    fig.savefig(out / f"{name}.png", dpi=155)
    plt.close(fig)


def motion_panel(out):
    with TemporaryDirectory() as temp:
        folder = Path(temp)
        files = []
        for i, scale in enumerate((0.012, 0.034)):
            steps = RNG.normal(0, scale, (250, 6))
            steps[:, :3] /= 50
            steps[[48, 100, 172], 3] += (0.18 + i * 0.18)
            path = folder / f"arbitrary-run-{i + 1}.par"
            np.savetxt(path, np.cumsum(steps, axis=0))
            files.append(path)
        motion(SimpleNamespace(supine=files[0], upright=files[1], out=folder))
        fig, axes = plt.subplots(2, 1, figsize=(8, 4.5), sharex=True, sharey=True,
                                 layout="constrained")
        for ax, posture, colour in zip(axes, ("supine", "upright"), (BLUE, ORANGE)):
            table = np.genfromtxt(folder / f"{posture}_fd.tsv", delimiter="\t", names=True)
            ax.plot(table["minutes_from_start"], table["power_fd_mm"], color=colour, linewidth=1.2)
            ax.axhline(0.2, color="#59636d", linestyle="--")
            ax.set_ylabel("Power FD (mm)")
            ax.set_title(f"Artificial {posture} motion", loc="left", fontsize=10)
    axes[-1].set_xlabel("Minutes from scan start")
    save(fig, out, "motion-example", "Motion QC with one scale and a 0.2 mm reference")


def spatial_panel(out):
    y, x = np.mgrid[-1:1:160j, -1:1:160j]
    brain = x**2 / 0.72**2 + y**2 / 0.92**2 <= 1
    wm = x**2 / 0.49**2 + y**2 / 0.65**2 <= 1
    grey = np.where(brain, 0.53 + 0.11 * np.cos(13 * x) * np.cos(12 * y), 0)
    anatomy = grey + wm * 0.30
    outside = (x**2 / 0.79**2 + y**2 / 0.96**2 <= 1) & ~brain
    anatomy[outside] = 0.18
    yy, xx = np.indices(anatomy.shape)
    epi = map_coordinates(grey + wm * 0.11, [yy, xx + 9 * np.exp(-((y + 0.42) / 0.28)**2)], order=1)
    epi[y < -0.68] = 0
    epi[(x**2 + (y + 0.50)**2 < 0.19**2)] *= 0.16
    usable = brain & (epi >= 0.25 * np.median(epi[brain]))
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 4.2), layout="constrained")
    for ax, image, title in zip(axes, (anatomy, epi, epi),
                                ("T1 + BET mask", "Mean EPI + T1 WM boundary", "EPI + usable region")):
        ax.imshow(image, origin="lower", cmap="gray", vmin=0, vmax=1)
        ax.set(title=title, xticks=[], yticks=[])
    axes[0].contour(brain, levels=[0.5], colors=[CYAN], linewidths=1.3)
    axes[1].contour(wm, levels=[0.5], colors=["#ffc475"], linewidths=1.4)
    axes[2].contour(usable, levels=[0.5], colors=[CYAN], linewidths=1.3)
    cutoff = np.flatnonzero(y[:, 0] < -0.68)[-1] + 0.5
    axes[2].axhline(cutoff, color="#ffc475", linestyle="--", linewidth=1)
    save(fig, out, "spatial-example", "Check masks, alignment, and signal coverage separately")


def ica_panel(out):
    x, y = np.mgrid[-1:1:120j, -1:1:120j]
    brain = x**2 / 0.77**2 + y**2 / 0.89**2 <= 1
    band = np.where(brain, np.exp(-((x + 0.20) / 0.045)**2), 0)
    ring = np.where(brain & (x**2 / 0.70**2 + y**2 / 0.82**2 > 1), 1, 0)
    time = np.arange(250) * TR
    trace = (np.cos(time / 19) + 0.42 * RNG.normal(0, 1, time.size))
    drift = RNG.normal(0, 0.45, time.size)
    drift[:3] += (5.2, 4.7, 2.1)
    fig, axes = plt.subplots(2, 2, figsize=(8.7, 6.5), layout="constrained")
    for row, active, series, title in zip(
        axes, (band, ring), (trace, drift), ("Narrow spatial band", "Peripheral spatial pattern")
    ):
        row[0].imshow(np.where(brain, 0.46, 0), cmap="gray", vmin=0, vmax=1)
        overlay = np.ma.masked_less(active, 0.15)
        row[0].imshow(overlay, cmap="autumn", alpha=0.85, vmin=0, vmax=1)
        row[0].set(title=title, xticks=[], yticks=[])
        row[1].plot(time / 60, series, color=BLUE, linewidth=1)
        row[1].axvspan(0, 3 * TR / 60, color=ORANGE, alpha=0.24)
        row[1].set(xlabel="Minutes", ylabel="Illustrative amplitude")
    save(fig, out, "ica-example", "Compare spatial morphology and time course before removing ICs")


def coverage_panel(out):
    supine = np.array([0.89, 0.91, 0.72, 0.69, 0.76, 0.86, 0.53, 0.98])
    upright = np.array([0.81, 0.88, 0.66, 0.82, 0.78, 0.80, 0.48, 0.95])
    native_min = np.array([96, 74, 57, 98, 34, 80, 31, 132])
    kept = (supine >= 0.7) & (upright >= 0.7) & (native_min >= 50)
    fig, ax = plt.subplots(figsize=(5.7, 5.6), layout="constrained")
    ax.scatter(supine[~kept], upright[~kept], color=ORANGE, label="Excluded")
    ax.scatter(supine[kept], upright[kept], color=BLUE, label="Kept in both")
    for i, (s, u, n) in enumerate(zip(supine, upright, native_min), 1):
        ax.annotate(f"{i}: {n}", (s, u), xytext=(2, 8), textcoords="offset points", fontsize=8)
    ax.axvline(0.7, color="#59636d", linestyle="--")
    ax.axhline(0.7, color="#59636d", linestyle="--")
    ax.set(xlabel="Supine usable fraction", ylabel="Upright usable fraction",
           xlim=(0.45, 1.06), ylim=(0.43, 1.04), aspect="equal")
    ax.legend(loc="lower right", title="Labels show ID : min native count", fontsize=8)
    save(fig, out, "coverage-example", "A parcel must pass both criteria in both runs")


def fc_panel(out):
    wm = RNG.normal(0, 1, 247)
    latent = RNG.normal(0, 1, (247, 2))
    nodes = latent[:, np.arange(8) % 2] * 0.7 + RNG.normal(0, 0.45, (247, 8))
    nodes += wm[:, None] * np.r_[np.full(123, 1.0), np.full(124, 1.7)][:, None]
    centred = wm - wm.mean()
    residual = nodes - np.outer(centred, centred @ (nodes - nodes.mean(axis=0)) / (centred @ centred))
    fig, axes = plt.subplots(2, 2, figsize=(7, 6.4), layout="constrained")
    for row, design, data in zip(axes, ("Before WM", "After WM"), (nodes, residual)):
        for ax, description, half in zip(row, ("first half", "second half"), HALVES):
            corr = np.corrcoef(data[half], rowvar=False)
            np.fill_diagonal(corr, np.nan)
            view = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
            ax.set(title=f"{design} · {description}", xticks=[], yticks=[])
    fig.colorbar(view, ax=axes, shrink=0.76, label="Pearson r")
    save(fig, out, "fc-example", "Keep node order and the Pearson scale fixed")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("figures"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for create in (motion_panel, spatial_panel, ica_panel, coverage_panel, fc_panel):
        create(args.out)


if __name__ == "__main__":
    main()
