#!/usr/bin/env bash
# run_uci_pipeline.sh — full UCI Cardiac CTP pipeline for one phase-2 rate.
#
# Steps per rate:
#   1. extract_peak_iodine_bae       (Bae AIF + per-tree contrast → .f32 + peak_metadata)
#   2. apply_contrast_at_peak        (V4 cross-product voxelization of trees into XCAT)
#   3. add_chambers_to_phantom       (patch chambers 19-22/24/25/28 with Bae C)
#   4. run_cta_sim_param             (BasisSimulator scan @ 120 kVp → FBP + HIR + DICOM)
#
# Usage:  bash run_uci_pipeline.sh  PHASE2_RATE  [SCAN_DELAY_S=6.1]
#   bash run_uci_pipeline.sh 2.5
#   bash run_uci_pipeline.sh 1.0
#
# Outputs go to:  phantom_ct_input/bae_uci_r{RATE}/{peak_iodine,peak_phantom,chamber_patched,cta_out}/
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0  PHASE2_RATE_ML_S  [SCAN_DELAY_S=6.1]"
    exit 1
fi
RATE="$1"
SCAN_DELAY="${2:-6.1}"

# Top-level paths (absolute)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.."  && pwd)"
TREE_DIR="${ROOT}/VascularTreeSim.jl/output"
PHANTOM_SRC="/home/molloi-lab/smb_mount/shared_drive/Shu Nie/PVAT_Analysis/digital phantoms/vmale50_1600x1400x500_8bit_little_endian_act_1.raw"

# Per-rate working dirs
WORK="${ROOT}/phantom_ct_input/bae_uci_r${RATE}"
PEAK_DIR="${WORK}/peak_iodine"
PHANTOM_PEAK_DIR="${WORK}/peak_phantom"
CHAMBER_DIR="${WORK}/chamber_patched"
CTA_OUT="${WORK}/cta_out"
mkdir -p "${PEAK_DIR}" "${PHANTOM_PEAK_DIR}" "${CHAMBER_DIR}" "${CTA_OUT}"

cd "${ROOT}"
echo "================================================================"
echo "UCI Cardiac CTP pipeline — phase-2 rate = ${RATE} mL/s, scan_delay = ${SCAN_DELAY} s"
echo "Working dir: ${WORK}"
echo "================================================================"

# ── Step 1: extract per-tree peak iodine with Bae AIF ────────────────
echo
echo "[1/4] extract_peak_iodine_bae"
echo "------------------------------------------------"
julia --project=FlowContrastSim.jl \
      FlowContrastSim.jl/scripts/extract_peak_iodine_bae.jl \
      "${TREE_DIR}" "${PEAK_DIR}" "${RATE}" "${SCAN_DELAY}"

# ── Step 2: V4 cross-product voxelization of trees into XCAT ─────────
echo
echo "[2/4] apply_contrast_at_peak (V4 cross-product, ~7 min on 64-core)"
echo "------------------------------------------------"
julia --project=VascularTreeSim.jl --threads=auto \
      VascularTreeSim.jl/scripts/apply_contrast_at_peak.jl \
      "${TREE_DIR}" "${PHANTOM_SRC}" "${PEAK_DIR}" "${PHANTOM_PEAK_DIR}"

# ── Step 3: chamber patching (XCAT 19/20/21/22/24/25/28 ← Bae C) ─────
echo
echo "[3/4] add_chambers_to_phantom"
echo "------------------------------------------------"
julia --project=VascularTreeSim.jl --threads=auto \
      VascularTreeSim.jl/scripts/add_chambers_to_phantom.jl \
      "${PHANTOM_PEAK_DIR}" "${PEAK_DIR}" "${CHAMBER_DIR}"

# ── Step 4: BasisSimulator scan @ 120 kVp ─────────────────────────────
echo
echo "[4/4] run_cta_sim_param @ 120 kVp"
echo "------------------------------------------------"
julia --project=phantom_ct_input/run_basis_sim \
      phantom_ct_input/run_basis_sim/run_cta_sim_param.jl \
      "${CHAMBER_DIR}" "${CTA_OUT}" --kvp 120 --mA 250

echo
echo "================================================================"
echo "Pipeline DONE for rate=${RATE} mL/s"
echo "  peak_iodine        : ${PEAK_DIR}"
echo "  peak_phantom       : ${PHANTOM_PEAK_DIR}"
echo "  chamber_patched    : ${CHAMBER_DIR}"
echo "  CT output (DICOM)  : ${CTA_OUT}"
echo "================================================================"
