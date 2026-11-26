"""Exploratory four-scan ROI connectivity and censoring sensitivity analysis."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy.ndimage import binary_erosion

from reproduce import atlas_lookup, atlas_names, deep_white_matter, nearest, on_grid, values

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True, help="Archive with the volunteer1 and volunteer2 directories")
parser.add_argument("--fsldir", type=Path, required=True, help="FSL installation containing the atlas")
parser.add_argument("--out", type=Path, required=True, help="New directory for QC and results")
args = parser.parse_args()
ROOT = args.root
FSL = args.fsldir / "data"
OUT = args.out
OUT.mkdir(parents=True)

def save_csv(name, header, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


atlas_img = nib.load(FSL / "atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz")
atlas_mni = np.asarray(atlas_img.dataobj)
names = atlas_names(FSL / "atlases/HarvardOxford-Cortical.xml")
label_ids = range(1, 49)
if set(names) != set(label_ids):
    raise ValueError("Expected the 48-label Harvard–Oxford cortical atlas")
config = {
    "v1_supine": {
        "base": "volunteer1/supine",
        "mean": "registration_qc/BOLD_mc_mean.nii.gz",
        "mask": "ica_qc/trim3/supine_trim3.ica/mask.nii.gz",
        "atlas": "connectivity_qc/brain_fit/HO_cort_BOLD.nii.gz",
        "mat": "connectivity_qc/brain_fit/MNI_to_BOLD_affine.mat",
        "motion": "motion_qc/BOLD_mc.par",
        "mc": "motion_qc/BOLD_mc.nii.gz",
        "hp": "ica_qc/trim3/BOLD_hp100.nii.gz",
        "clean": "ica_qc/trim3/slice_denoise_qc/BOLD_hp100_regfilt_ic1_2_3.nii.gz",
        "wm": "connectivity_qc/tissue_qc/WM_pve_in_BOLD.nii.gz",
    },
    "v1_upright": {
        "base": "volunteer1/upright",
        "mean": "func/BOLD_mean.nii.gz",
        "mask": "func/upright.ica/mask.nii.gz",
        "atlas": "reg/atlas_BOLD.nii.gz",
        "mat": "reg/MNI_to_BOLD.mat",
        "motion": "func/BOLD_mc.par",
        "mc": "func/BOLD_mc.nii.gz",
        "hp": "func/BOLD_hp100.nii.gz",
        "clean": "func/BOLD_ic_clean.nii.gz",
        "wm": "reg/WM_BOLD.nii.gz",
    },
    "v2_supine": {
        "base": "volunteer2/supine",
        "mean": "analysis_qc/BOLD_mean_mc.nii.gz",
        "mask": "analysis_qc/single.ica/mask.nii.gz",
        "atlas": "analysis_qc/atlas_native.nii.gz",
        "mat": "coverage_qc/MNI_to_BOLD_centered.mat",
        "motion": "analysis_qc/BOLD_mc.par",
        "mc": "analysis_qc/BOLD_trim3.nii.gz",
        "hp": "analysis_qc/BOLD_hp100.nii.gz",
        "clean": "analysis_qc/BOLD_ica_clean.nii.gz",
        "wm": "analysis_qc/WM_pve_native.nii.gz",
    },
    "v2_upright": {
        "base": "volunteer2/upright",
        "mean": "analysis_qc/BOLD_mean_mc.nii.gz",
        "mask": "analysis_qc/single.ica/mask.nii.gz",
        "atlas": "analysis_qc/atlas_native.nii.gz",
        "mat": "coverage_qc/MNI_to_BOLD.mat",
        "motion": "analysis_qc/BOLD_mc.par",
        "mc": "analysis_qc/BOLD_trim3.nii.gz",
        "hp": "analysis_qc/BOLD_hp100.nii.gz",
        "clean": "analysis_qc/BOLD_ica_clean.nii.gz",
        "wm": "analysis_qc/WM_binary_native.nii.gz",
    },
}
runs = {}
coverage_rows = []
for name, cfg in config.items():
    folder = ROOT / cfg["base"]
    ref = nib.load(folder / cfg["mean"])
    mean = values(ref)
    brain = on_grid(folder / cfg["mask"], ref) > 0
    labels = on_grid(folder / cfg["atlas"], ref).astype(int)
    if not brain.any() or not set(np.unique(labels)) <= {0, *label_ids}:
        raise ValueError(f"Invalid mask or atlas labels: {name}")
    usable = brain & (mean >= .25 * np.median(mean[brain]))
    transform = atlas_lookup(ref, atlas_img, np.loadtxt(folder / cfg["mat"]))
    in_atlas = nearest(usable, transform, atlas_mni.shape)
    coverage = np.array([in_atlas[atlas_mni == i].mean() for i in range(1,49)])
    counts = np.array([np.count_nonzero(usable & (labels == i)) for i in range(1,49)])
    for i in range(48):
        coverage_rows.append([name, i + 1, names[i + 1], coverage[i], counts[i], coverage[i] >= .7 and counts[i] >= 50])
    motion = np.loadtxt(folder / cfg["motion"])
    mc = on_grid(folder / cfg["mc"], ref)
    if motion.shape[0] == 250:
        motion = motion[3:]
    if mc.shape[-1] == 250:
        mc = mc[..., 3:]
    if motion.shape != (247, 6) or mc.shape != (*mean.shape, 247):
        raise ValueError(f"Expected 247 volumes and matching motion parameters: {name}")
    diff = np.abs(np.diff(motion, axis=0))
    # The first analysed volume has no preceding volume in the trimmed series.
    fd = np.r_[0., diff[:, 3:].sum(axis=1) + 50 * diff[:, :3].sum(axis=1)]
    erosion = binary_erosion(usable)
    sig = mc[erosion]
    if not len(sig) or not np.isfinite(sig).all() or np.median(sig) <= 0:
        raise ValueError(f"Invalid eroded signal mask for DVARS: {name}")
    dv = np.r_[0., np.sqrt(np.mean(np.diff(sig, axis=1) ** 2, axis=0)) / np.median(sig) * 100]
    limit = np.quantile(dv[1:], .75) + 1.5 * (np.quantile(dv[1:], .75) - np.quantile(dv[1:], .25))
    wm_prob = on_grid(folder / cfg["wm"], ref)
    wm_mask = deep_white_matter(wm_prob, usable, ref.header.get_zooms()[:3])
    wm_mask &= labels == 0
    print(name, "FD mean", fd[1:].mean(), "DVARS fence", limit, "WM", wm_mask.sum(), "coverage", sum((coverage >= .7) & (counts >= 50)), flush=True)
    if wm_mask.sum() < 50:
        raise ValueError(f"Insufficient deep WM voxels in {name}")
    runs[name] = dict(folder=folder, cfg=cfg, reference=ref, usable=usable, labels=labels,
                      coverage=coverage, counts=counts, fd=fd, dv=dv, dvlim=limit,
                      motion=motion, wm=wm_mask)
    save_csv(name + "_motion.csv", ["original_volume", "FD_mm", "DVARS_percent", "DVARS_fence"], zip(range(4,251),fd,dv,np.repeat(limit,247)))
    del mc, sig

save_csv("coverage_four_scans.csv", ["scan", "label", "region", "coverage", "usable_voxels", "passes"], coverage_rows)
shared = [i for i in range(1,49) if all(r["coverage"][i-1] >= .7 and r["counts"][i-1] >= 50 for r in runs.values())]
print("Common nodes", len(shared), shared, flush=True)
if len(shared) < 2:
    raise ValueError("Fewer than two common nodes")
save_csv("common_nodes.csv", ["label", "region"], [(i,names[i]) for i in shared])

for name, r in runs.items():
    r["series"] = {}
    for key in ("hp", "clean"):
        data = on_grid(r["folder"] / r["cfg"][key], r["reference"])
        if data.shape != (*r["usable"].shape, 247):
            raise ValueError(f"Expected a 247-volume {key} derivative: {name}")
        r["series"][key] = np.column_stack([data[r["usable"] & (r["labels"] == i)].mean(axis=0) for i in shared])
        r["series"][key + "_wm"] = data[r["wm"]].mean(axis=0)
        del data

conditions = ("hp_all", "ica_all", "motion_wm_all", "motion_wm_fd03", "motion_wm_fd03_dvars", "motion_wm_matched")
retained = {k: (r["fd"] <= .3) & (r["dv"] <= r["dvlim"]) for k,r in runs.items()}
match_n = min(int(v.sum()) for v in retained.values())
print("Matched retained volumes",match_n,flush=True)
if match_n < 30:
    print("WARNING fewer than 30 volumes in matched analysis",flush=True)
matrices = {}
metrics = []
split_rows = []
for key,r in runs.items():
    X = np.column_stack([np.ones(247), r["motion"], np.vstack([np.zeros(6),np.diff(r["motion"],axis=0)]), r["series"]["clean_wm"]])
    X = np.nan_to_num(X)
    selections = {"hp_all":np.ones(247,dtype=bool), "ica_all":np.ones(247,dtype=bool), "motion_wm_all":np.ones(247,dtype=bool),
                  "motion_wm_fd03":r["fd"]<=.3, "motion_wm_fd03_dvars":retained[key]}
    ids=np.flatnonzero(retained[key])
    chosen=ids[np.linspace(0,len(ids)-1,match_n,dtype=int)]
    selections["motion_wm_matched"]=np.isin(np.arange(247),chosen)
    for condition, keep in selections.items():
        y=r["series"]["hp"] if condition == "hp_all" else r["series"]["clean"]
        if condition not in ("hp_all", "ica_all"):
            beta=np.linalg.lstsq(X[keep],y[keep],rcond=None)[0]
            y=y-X@beta
        mat=np.corrcoef(y[keep],rowvar=False)
        if not np.isfinite(mat).all():
            raise ValueError(f"Non-finite matrix: {key} {condition}")
        matrices[(key,condition)]=mat
        keep_ids = np.flatnonzero(keep)
        left = keep_ids[keep_ids < 123]
        right = keep_ids[keep_ids >= 124]
        if len(left) >= 20 and len(right) >= 20:
            first = np.corrcoef(y[left], rowvar=False)
            second = np.corrcoef(y[right], rowvar=False)
            tri = np.triu_indices(len(shared), 1)
            agreement = np.corrcoef(first[tri], second[tri])[0, 1]
            split_rows.append([key, condition, len(left), len(right), agreement])
        else:
            split_rows.append([key, condition, len(left), len(right), ""])
        np.savetxt(OUT / (key+"_"+condition+".csv"),mat,delimiter=",",fmt="%.6f",
                   header=",".join(str(i) for i in shared),comments="")
        metrics.append([key,condition,int(keep.sum()),np.median(mat[np.triu_indices(len(shared),1)]),np.linalg.matrix_rank(X[keep]) if condition not in ("hp_all", "ica_all") else 0])
    save_csv(key+"_retained.csv",["original_volume","fd03_dvars","matched"],zip(range(4,251),retained[key].astype(int),selections["motion_wm_matched"].astype(int)))
save_csv("analysis_summary.csv",["scan","condition","retained","median_edge_r","nuisance_design_rank"],metrics)
save_csv("split_half_sensitivity.csv",["scan","condition","first_count","second_count","edge_pattern_correlation"],split_rows)

edge=np.triu_indices(len(shared),1)
all_changes = []
for condition in conditions:
    effects=[]
    for v in ("v1","v2"):
        effects.append(matrices[(v+"_upright",condition)]-matrices[(v+"_supine",condition)])
    np.savez_compressed(OUT/(condition+"_matrices.npz"),v1_diff=effects[0],v2_diff=effects[1],labels=shared)
    fig,axes=plt.subplots(1,2,figsize=(12,5.5),constrained_layout=True)
    for ax,diff,v in zip(axes,effects,("Volunteer 1","Volunteer 2")):
        plot=ax.imshow(diff,cmap="RdBu_r",vmin=-1,vmax=1)
        ax.set_title(v+"  upright minus supine")
        ax.set_xticks(range(len(shared)),shared,rotation=90,fontsize=6)
        ax.set_yticks(range(len(shared)),shared,fontsize=6)
    fig.colorbar(plot,ax=axes,shrink=.7,label="Difference in Pearson r")
    fig.suptitle(condition+" | identical atlas label order")
    fig.savefig(OUT/(condition+"_difference.png"),dpi=150)
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,10),constrained_layout=True)
    for ax,scan in zip(axes.flat,("v1_supine","v1_upright","v2_supine","v2_upright")):
        im=ax.imshow(matrices[(scan,condition)],cmap="RdBu_r",vmin=-1,vmax=1)
        ax.set_title(scan)
        ax.set_xticks(range(len(shared)),shared,rotation=90,fontsize=5)
        ax.set_yticks(range(len(shared)),shared,fontsize=5)
    fig.colorbar(im,ax=axes,shrink=.7,label="Pearson r")
    fig.suptitle(condition+" | shared atlas labels")
    fig.savefig(OUT/(condition+"_four_matrices.png"),dpi=150)
    plt.close(fig)
    a,b=(e[edge] for e in effects)
    all_changes.append(np.stack([a, b]))
    same=np.sign(a)==np.sign(b)
    save_csv(condition+"_edge_differences.csv",["label1","region1","label2","region2","v1_difference","v2_difference","same_direction"],
             [(shared[i],names[shared[i]],shared[j],names[shared[j]],a[k],b[k],bool(same[k])) for k,(i,j) in enumerate(zip(*edge))])
    print(condition,"edge effects same sign",same.sum(),"/",len(same),"effect correlation",np.corrcoef(a,b)[0,1],flush=True)

changes = np.stack(all_changes)
directions = np.sign(changes)
consistent = np.all(directions == directions[0, 0], axis=(0, 1))
consistent_rows = []
for index in np.flatnonzero(consistent):
    i, j = edge[0][index], edge[1][index]
    a, b = shared[i], shared[j]
    magnitudes = np.abs(changes[:, :, index])
    consistent_rows.append((a, names[a], b, names[b],
                            "increase" if directions[0, 0, index] > 0 else "decrease",
                            f"{magnitudes.min():.6f}", f"{magnitudes.max():.6f}"))
save_csv("edges_consistent_all_settings.csv",
         ["label1", "region1", "label2", "region2", "direction",
          "minimum_absolute_change", "maximum_absolute_change"], consistent_rows)
print("Same direction across all six conditions in both participants:", sum(consistent), flush=True)

fig,axes=plt.subplots(4,2,figsize=(13,9),sharex=True,constrained_layout=True)
for row,(key,r) in enumerate(runs.items()):
    axes[row,0].plot(np.arange(4,251),r["fd"],lw=.8)
    axes[row,0].axhline(.3,color="darkorange",ls="--")
    axes[row,0].set_ylabel(key+"\nFD mm")
    axes[row,1].plot(np.arange(4,251),r["dv"],lw=.8,color="teal")
    axes[row,1].axhline(r["dvlim"],color="darkorange",ls="--")
    axes[row,1].set_ylabel("DVARS %")
axes[0,0].set_title("FD after first three volumes")
axes[0,1].set_title("DVARS after motion correction")
axes[-1,0].set_xlabel("Original volume")
axes[-1,1].set_xlabel("Original volume")
fig.savefig(OUT/"four_scan_motion_dvars.png",dpi=150)
plt.close(fig)
print("Results",OUT,flush=True)
