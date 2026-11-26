#!/usr/bin/env bash
# Reproduce the Volunteer 1 pilot in separate, inspectable stages.
set -euo pipefail

usage() {
  cat <<'HELP'
Usage: bash run.sh STAGE WORK [arguments]

  check              Validate input headers and record checksums and versions.
  motion             Create raw means, MCFLIRT outputs, BOLD masks and FD.
  register           Run BBR and map the atlas and white matter to native BOLD.
  ica                Run the 250-volume reference and both trimmed ICA runs.
  denoise S_IC U_IC   Regress three reviewed IDs per posture, e.g. 1,2,3 2,3,14.
  analyse            Recompute ICA checks, slice QC and split-half connectivity.
  figures            Draw teaching overlays and plots from the existing outputs.
  replay ARCHIVE     Link Volunteer 1 archived inputs in a new WORK and run QC.
  coverage-replay ARCHIVE PARCELS_CSV
                     Link Volunteer 2 archived inputs and recompute coverage.

Configure FSLDIR and set PILOT_PYTHON to an interpreter with NumPy, SciPy,
NiBabel and Matplotlib. WORK is a fresh output directory for the replay stages.
HELP
}

die() { printf 'FSL pilot: %s\n' "$*" >&2; exit 1; }
require() {
  for path in "$@"; do
    [[ -r "$path" ]] || die "Missing input: $path"
  done
}

[[ $# -gt 0 && "$1" != --help && "$1" != -h ]] || { usage; exit 0; }
[[ $# -ge 2 ]] || die 'Specify a stage and work directory. See --help.'
stage=$1
work=$2
shift 2
case "$stage" in
  check|motion|register|ica|analyse|figures) expected=0 ;;
  replay) expected=1 ;;
  denoise|coverage-replay) expected=2 ;;
  *) die "Unknown stage: $stage" ;;
esac
[[ $# -eq "$expected" ]] || die "Expected $expected additional arguments for $stage."
if [[ "$stage" == replay || "$stage" == coverage-replay ]]; then
  [[ ! -e "$work" && ! -L "$work" ]] || die "Choose a new replay directory: $work"
fi

script_dir=$(cd "$(dirname "$0")" && pwd -P)
export FSLDIR="${FSLDIR:-$HOME/fsl}"
require "$FSLDIR/etc/fslconf/fsl.sh"
# FSL's initialisation may read unset variables.
set +u
source "$FSLDIR/etc/fslconf/fsl.sh"
set -u
export PATH="$FSLDIR/share/fsl/bin:$PATH" FSLOUTPUTTYPE=NIFTI_GZ
export OMP_NUM_THREADS="${PILOT_THREADS:-2}" OPENBLAS_NUM_THREADS="${PILOT_THREADS:-2}"
export MKL_NUM_THREADS="$OPENBLAS_NUM_THREADS"
python="${PILOT_PYTHON:-python3}"
command -v "$python" >/dev/null || die "Python not found: $python"
script="$script_dir/reproduce.py"
standard="$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"
atlas="$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz"
labels="$FSLDIR/data/atlases/HarvardOxford-Cortical.xml"

check_inputs() {
  require "$work"/{supine,upright}/{t1,bold_raw}.nii.gz
  for posture in supine upright; do
    t1="$work/$posture/t1.nii.gz"
    bold="$work/$posture/bold_raw.nii.gz"
    [[ $(fslval "$t1" dim4) -eq 1 ]] || die "$posture T1 should be 3D."
    [[ $(fslval "$bold" dim4) -eq 250 ]] || die "$posture BOLD should have 250 volumes."
    fslval "$bold" pixdim4 | awk '{d=$1-3.8; exit !(d*d < 0.000001)}' ||
      die "$posture BOLD should have TR 3.8 seconds. Check the time units."
    fslinfo "$bold"
  done
  sha256sum "$work"/{supine,upright}/{t1,bold_raw}.nii.gz > "$work/qc/input_sha256.txt"
}

motion() {
  check_inputs
  for posture in supine upright; do
    (
      cd "$work/$posture"
      fslmaths bold_raw -Tmean raw_mean
      mcflirt -in bold_raw -out mc -plots -rmsrel -rmsabs
      fslmaths mc -Tmean mean
      bet mean epi_brain -f 0.30 -m
    )
  done
  "$python" "$script" motion --supine "$work/supine/mc.par" \
    --upright "$work/upright/mc.par" --out "$work/qc"
}

register() {
  require "$standard" "$atlas" "$work"/{supine,upright}/{t1,mean}.nii.gz
  local -a affine_fit=(-dof 12 -cost corratio -searchrx -25 25 -searchry -25 25 -searchrz -25 25)
  (
    cd "$work/supine"
    bet t1 t1_brain -f 0.37 -m
    epi_reg --epi=mean --t1=t1 --t1brain=t1_brain --out=bold_to_t1 --noclean
    flirt -in t1_brain -ref t1_brain -applyisoxfm 1 -out t1_brain_1mm
    flirt -in t1_brain_1mm -ref "$standard" -omat iso_to_mni.mat -out t1_iso_mni "${affine_fit[@]}"
    "$python" "$script" bridge --native t1.nii.gz --isotropic t1_brain_1mm.nii.gz \
      --iso-to-mni iso_to_mni.mat --output t1_to_mni.mat
  )
  (
    cd "$work/upright"
    flirt -in t1 -ref t1 -applyisoxfm 1 -out t1_1mm
    bet t1_1mm t1_brain -R -f 0.37 -m
    fast -o t1_fast t1_brain
    fslmaths t1_fast_pve_2 -thr 0.5 -bin t1_wm
    epi_reg --epi=mean --t1=t1_1mm --t1brain=t1_brain --wmseg=t1_wm --out=bold_to_t1
    flirt -in t1_brain -ref "$standard" -omat t1_to_mni.mat -out t1_iso_mni "${affine_fit[@]}"
  )
  for posture in supine upright; do
    (
      cd "$work/$posture"
      convert_xfm -omat bold_to_mni.mat -concat t1_to_mni.mat bold_to_t1.mat
      convert_xfm -omat mni_to_bold.mat -inverse bold_to_mni.mat
      convert_xfm -omat t1_to_bold.mat -inverse bold_to_t1.mat
      flirt -in "$atlas" -ref mean -applyxfm -init mni_to_bold.mat \
        -interp nearestneighbour -out atlas_native
      if [[ "$posture" == supine ]]; then wm=bold_to_t1_fast_pve_2; else wm=t1_fast_pve_2; fi
      flirt -in "$wm" -ref mean -applyxfm -init t1_to_bold.mat \
        -interp trilinear -out wm_native
    )
  done
}

ica() {
  require "$work"/{supine,upright}/{mc,mean,epi_brain_mask}.nii.gz
  for directory in "$work"/supine/{initial,trim3}.ica "$work/upright/trim3.ica"; do
    [[ ! -e "$directory" && ! -L "$directory" ]] || die "ICA output already exists: $directory"
  done
  export OMP_NUM_THREADS="${PILOT_ICA_THREADS:-4}" OPENBLAS_NUM_THREADS="${PILOT_ICA_THREADS:-4}"
  export MKL_NUM_THREADS="$OPENBLAS_NUM_THREADS"
  (
    cd "$work/supine"
    fslmaths mc -bptf 13.157895 -1 -add mean hp100_all -odt float
    melodic -i hp100_all.nii.gz -o initial.ica --tr=3.8 --nobet \
      -m epi_brain_mask.nii.gz --bgimage=mean.nii.gz --report --Ostats --Omean
  )
  for posture in supine upright; do
    (
      cd "$work/$posture"
      fslroi mc trim3 3 247
      fslmaths trim3 -Tmean trim3_mean
      fslmaths trim3 -bptf 13.157895 -1 -add trim3_mean hp100 -odt float
      melodic -i hp100.nii.gz -o trim3.ica --tr=3.8 --nobet \
        -m epi_brain_mask.nii.gz --Omean
    )
  done
}

validate_ids() {
  local ids=$1 mix=$2
  [[ "$ids" =~ ^[1-9][0-9]*,[1-9][0-9]*,[1-9][0-9]*$ ]] ||
    die "Expected three comma-separated 1-based IC IDs: $ids"
  awk -v ids="$ids" '
    NR == 1 {
      split(ids, selected, ",")
      for (i in selected) {
        id = selected[i] + 0
        if (id > NF || seen[id]++) bad = 1
      }
    }
    END { exit (NR != 247 || bad) }
  ' "$mix" || die "Check distinct IDs, ICA columns, and 247 rows in $mix."
}

denoise() {
  require "$work"/{supine,upright}/{hp100.nii.gz,trim3.ica/melodic_mix,trim3.ica/mask.nii.gz}
  validate_ids "$1" "$work/supine/trim3.ica/melodic_mix"
  validate_ids "$2" "$work/upright/trim3.ica/melodic_mix"
  printf 'posture\tselected_ics\nsupine\t%s\nupright\t%s\n' "$1" "$2" > "$work/qc/selected_ics.tsv"
  for posture in supine upright; do
    if [[ "$posture" == supine ]]; then ids=$1; else ids=$2; fi
    (
      cd "$work/$posture"
      fsl_regfilt -i hp100.nii.gz -d trim3.ica/melodic_mix \
        -m trim3.ica/mask.nii.gz -f "$ids" -o clean.nii.gz
      [[ $(fslval clean dim4) -eq 247 ]] || die "$posture clean BOLD should have 247 volumes."
    )
  done
}

analyse() {
  require "$atlas" "$labels" "$work"/{supine,upright}/{mc.par,hp100.nii.gz,clean.nii.gz}
  "$python" "$script" motion --supine "$work/supine/mc.par" \
    --upright "$work/upright/mc.par" --out "$work/qc"
  "$python" "$script" ica --initial "$work/supine/initial.ica" \
    --trimmed "$work/supine/trim3.ica" --out "$work/qc"
  "$python" "$script" slice-qc --before "$work/supine/hp100.nii.gz" \
    --after "$work/supine/clean.nii.gz" --mask "$work/supine/trim3.ica/mask.nii.gz" \
    --out "$work/qc"
  "$python" "$script" connectivity --supine "$work/supine" --upright "$work/upright" \
    --atlas "$atlas" --labels "$labels" --out "$work/qc"
}

figures() {
  "$python" "$script_dir/plot_qc.py" "$work"
}

link_input() {
  require "$1"
  ln -s -- "$1" "$2"
}

replay() {
  local source_dir s u src_path out_path
  source_dir=$(cd "$1" && pwd -P)
  s="$source_dir/supine"
  u="$source_dir/upright"
  mkdir "$work/supine" "$work/upright"
  while read -r src_path out_path; do
    link_input "$s/$src_path" "$work/supine/$out_path"
  done <<'SUPINE'
BOLD_original.nii.gz bold_raw.nii.gz
T1_original.nii.gz t1.nii.gz
motion_qc/BOLD_mc.par mc.par
registration_qc/BOLD_mc_mean.nii.gz mean.nii.gz
registration_qc/T1_brain_f037.nii.gz t1_brain.nii.gz
registration_qc/T1_brain_f037_mask.nii.gz t1_brain_mask.nii.gz
ica_qc/BOLD_bet_mask.nii.gz epi_brain_mask.nii.gz
ica_qc/supine_single.ica initial.ica
ica_qc/trim3/supine_trim3.ica trim3.ica
ica_qc/trim3/BOLD_hp100.nii.gz hp100.nii.gz
ica_qc/trim3/slice_denoise_qc/BOLD_hp100_regfilt_ic1_2_3.nii.gz clean.nii.gz
connectivity_qc/brain_fit/HO_cort_BOLD.nii.gz atlas_native.nii.gz
connectivity_qc/brain_fit/MNI_to_BOLD_affine.mat mni_to_bold.mat
connectivity_qc/tissue_qc/WM_pve_in_BOLD.nii.gz wm_native.nii.gz
SUPINE
  while read -r src_path out_path; do
    link_input "$u/$src_path" "$work/upright/$out_path"
  done <<'UPRIGHT'
BOLD_original.nii.gz bold_raw.nii.gz
T1_original.nii.gz t1.nii.gz
func/BOLD_mc.par mc.par
func/BOLD_mean.nii.gz mean.nii.gz
reg/T1_brain.nii.gz t1_brain.nii.gz
reg/T1_brain_mask.nii.gz t1_brain_mask.nii.gz
func/BOLD_bet_mask.nii.gz epi_brain_mask.nii.gz
func/upright.ica trim3.ica
func/BOLD_hp100.nii.gz hp100.nii.gz
func/BOLD_ic_clean.nii.gz clean.nii.gz
reg/atlas_BOLD.nii.gz atlas_native.nii.gz
reg/MNI_to_BOLD.mat mni_to_bold.mat
reg/WM_BOLD.nii.gz wm_native.nii.gz
UPRIGHT
  check_inputs
  printf 'posture\tselected_ics\nsupine\t1,2,3\nupright\t2,3,14\n' > "$work/qc/selected_ics.tsv"
  analyse
  figures
}

coverage_replay() {
  local source_dir src mat
  source_dir=$(cd "$1" && pwd -P)
  require "$2" "$atlas"
  mkdir "$work/supine" "$work/upright"
  for posture in supine upright; do
    src="$source_dir/$posture/coverage_qc"
    if [[ "$posture" == supine ]]; then mat=MNI_to_BOLD_centered.mat; else mat=MNI_to_BOLD.mat; fi
    link_input "$src/BOLD_mean.nii.gz" "$work/$posture/mean.nii.gz"
    link_input "$src/BOLD_bet_mask.nii.gz" "$work/$posture/epi_brain_mask.nii.gz"
    link_input "$src/$mat" "$work/$posture/mni_to_bold.mat"
  done
  "$python" "$script" coverage --supine "$work/supine" --upright "$work/upright" \
    --atlas "$atlas" --pilot-parcels "$2" --out "$work/four_scan_coverage.csv"
}

mkdir -p -- "$work"
work=$(cd "$work" && pwd -P)
case "$stage" in
  replay|coverage-replay) touch "$work/.archived-input-links" ;;
  motion|register|ica|denoise)
    [[ ! -e "$work/.archived-input-links" ]] || die 'Use a separate workspace for preprocessing.'
    ;;
esac
mkdir -p "$work/qc" "$work/logs"
log=$(mktemp "$work/logs/$stage.XXXXXXXX")
{
  printf 'Started: %s\n' "$(date -Iseconds)"
  printf 'Code revision: %s\n' "$(git -C "$script_dir" rev-parse HEAD 2>/dev/null || printf unknown)"
  sha256sum "$script_dir"/{run.sh,reproduce.py,plot_qc.py}
  printf 'FSLDIR: %s\n' "$FSLDIR"
  if [[ -r "$FSLDIR/etc/fslversion" ]]; then cat "$FSLDIR/etc/fslversion"; fi
  for entry in "$FSLDIR"/conda-meta/fsl-{base,bet2,mcflirt,flirt,fast4,avwutils,melodic,data_standard,data_atlases}-*.json; do
    if [[ -f "$entry" ]]; then basename "$entry" .json; fi
  done
  "$python" -c 'import sys, numpy, scipy, nibabel, matplotlib; print(sys.version); print({m.__name__: m.__version__ for m in (numpy, scipy, nibabel, matplotlib)})'
  printf 'Command: '; printf '%q ' "$stage" "$work" "$@"; printf '\n'
  (
    set -x
    case "$stage" in
      check) check_inputs ;;
      coverage-replay) coverage_replay "$@" ;;
      *) "$stage" "$@" ;;
    esac
  )
  printf 'Completed: %s\n' "$(date -Iseconds)"
} 2>&1 | tee "$log"
