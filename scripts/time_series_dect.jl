#!/usr/bin/env julia
# time_series_dect.jl — sweep BasisSimulator 120 kVp recon every Δt seconds
# from injection start, sampling chamber HU at each time. Chamber labels
# are re-patched in-memory per time step from the Bae V2 PBPK chamber
# concentration at that t (with V2's variable-volumetric-flow boost).
#
# IMPORTANT design choice — extended chamber-only bins.
#   We reuse the V1 rate=1.5 100kg phantom which has tree labels 256..10355
#   built with iodine_max=7.93 mgI/mL. V2 chamber peaks (RV ~12 mgI/mL) blow
#   past that, so re-using the same bins would saturate. Instead we DEFINE
#   a NEW chamber-only label band 10356..10456 with iodine_max=13 mgI/mL
#   (101 bins, all 100%-blood). At each time step, chamber voxels (originally
#   patched at labels 10342..10350) are re-stamped to this new band. Tree
#   labels (256..10355) are untouched.
#
# Per-time iteration cost: ~50–60 s on RTX 4090 (forward project + FBP).
# Only the LAST time-step DICOM is saved. All sampled HU accumulate in
# chamber_hu_timeseries.csv. Plot is rendered separately from the CSV.
#
# CLI:
#   julia --project=. time_series_dect.jl PHANTOM_DIR OUTPUT_DIR \
#         [--rate 1.5] [--dt 1.0] [--t_end 60.0] [--weight 100]

import Pkg
# This file lives in FlowContrastSim.jl/scripts/ but USES BasisSimulator, so we
# activate the BS scripts project (which has both BasisSimulator and a Pkg.dev'd
# FlowContrastSim) instead of FCS itself.
const _BS_SCRIPTS = joinpath(@__DIR__, "..", "..", "phantom_ct_input", "run_basis_sim")
Pkg.activate(_BS_SCRIPTS)

import BasisSimulator as BS
import Unitful
using Unitful: @u_str
using Statistics: mean
using TOML
using Printf
using Dates: now
using FlowContrastSim   # Bae V2 (Q_central = Q_total + Q_inj boost)

include(joinpath(_BS_SCRIPTS, "write_dicom.jl"))

# ── 1. CLI ─────────────────────────────────────────────────────────────
function parse_args(args)
    pos = String[]
    rate = 1.5
    dt = 1.0
    t_end = 60.0
    weight_kg = 100.0
    height_cm = 173.0
    phase1_vol = 66.0;  phase1_rate = 5.0
    phase2_vol = 33.0
    phase3_vol = 30.0;  phase3_rate = 2.5
    conc = 370.0
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--rate"
            rate = parse(Float64, args[i+1]); i += 2
        elseif a == "--dt"
            dt = parse(Float64, args[i+1]); i += 2
        elseif a == "--t_end"
            t_end = parse(Float64, args[i+1]); i += 2
        elseif a == "--weight"
            weight_kg = parse(Float64, args[i+1]); i += 2
        else
            push!(pos, a); i += 1
        end
    end
    length(pos) >= 2 || error(
        "Usage: time_series_dect.jl PHANTOM_DIR OUTPUT_DIR [--rate RATE] [--dt DT] [--t_end T_END]")
    return (phantom_dir=abspath(pos[1]), output_dir=abspath(pos[2]),
            phase2_rate=rate, dt=dt, t_end=t_end,
            weight_kg=weight_kg, height_cm=height_cm, conc=conc,
            phase1_vol=phase1_vol, phase1_rate=phase1_rate,
            phase2_vol=phase2_vol,
            phase3_vol=phase3_vol, phase3_rate=phase3_rate)
end
const A = parse_args(ARGS)
isdir(A.phantom_dir) || error("phantom_dir not found: $(A.phantom_dir)")
mkpath(A.output_dir)
const KVP = 120.0
const MA = 250.0
const VIEWS = 500
const SCAN_TIMES = collect(0.0:A.dt:A.t_end)
@info "time_series_dect" phantom_dir=A.phantom_dir output_dir=A.output_dir phase2_rate=A.phase2_rate dt=A.dt t_end=A.t_end n_time_points=length(SCAN_TIMES)

# ── 2. Bae V2 once ─────────────────────────────────────────────────────
@info "Running Bae V2 PBPK"
patient = Patient(weight_kg=A.weight_kg, height_cm=A.height_cm)
protocol = TriphasicProtocol(
    weight_kg=A.weight_kg, contrast_concentration_mgI_ml=A.conc,
    phase1_volume_ml=A.phase1_vol, phase1_rate_ml_s=A.phase1_rate,
    phase2_volume_ml=A.phase2_vol, phase2_rate_ml_s=A.phase2_rate,
    phase3_volume_ml=A.phase3_vol, phase3_rate_ml_s=A.phase3_rate,
    phase2_dilution=1.0, phase3_dilution=0.0,
)
bae = simulate_central_circulation(patient, protocol; tspan=(0.0, A.t_end + 5.0), dt_save=0.05)
@printf("  CO = %.0f mL/min  (Bae formula for %g kg / %g cm)\n",
        FlowContrastSim.cardiac_output_ml_min(patient), A.weight_kg, A.height_cm)
for (name, C) in (("aorta",bae.C_aorta),("RV",bae.C_RV),("LV",bae.C_LV),
                  ("PA",bae.C_pulm_artery),("PV",bae.C_pulm_vein))
    pk = maximum(C); pt = bae.times[argmax(C)]
    @printf("  Bae V2 %-5s peak %.2f mgI/mL @ t=%.1fs\n", name, pk, pt)
end

# Dump Bae traces to CSV for record
open(joinpath(A.output_dir, "bae_v2_traces.csv"), "w") do io
    println(io, "t_s,C_aorta,C_RV,C_LV,C_pulm_artery,C_pulm_vein")
    for k in eachindex(bae.times)
        @printf(io, "%.3f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
            bae.times[k], bae.C_aorta[k], bae.C_RV[k], bae.C_LV[k],
            bae.C_pulm_artery[k], bae.C_pulm_vein[k])
    end
end

function c_at(times, C, t)
    t <= times[1] && return 0.0
    t >= times[end] && return C[end]
    i = searchsortedlast(times, t)
    frac = (t - times[i]) / (times[i+1] - times[i])
    return C[i]*(1-frac) + C[i+1]*frac
end

# ── 3. Phantom: load + downsample ─────────────────────────────────────
const VOXEL_SIZE_CM = (0.02, 0.02, 0.02)
const DOWNSAMPLE = 2

@info "Loading chamber-patched phantom"
cfg = TOML.parsefile(joinpath(A.phantom_dir, "phantom_manifest.toml"))
raw_path = joinpath(A.phantom_dir, cfg["phantom"]["raw_path"])
dims = Tuple(Int.(cfg["phantom"]["dims"]))
phantom_full = Array{UInt16}(undef, dims...)
t0 = time()
open(raw_path, "r") do io; read!(io, phantom_full); end
phantom_full = reverse(phantom_full; dims=(2,3))
nx_full, ny_full, nz_full = size(phantom_full)
nx_ds = nx_full ÷ DOWNSAMPLE; ny_ds = ny_full ÷ DOWNSAMPLE; nz_ds = nz_full ÷ DOWNSAMPLE
h = DOWNSAMPLE ÷ 2 + 1
phantom_labeled = Array{UInt16}(undef, nx_ds, ny_ds, nz_ds)
@inbounds for k in 1:nz_ds, j in 1:ny_ds, i in 1:nx_ds
    phantom_labeled[i,j,k] = phantom_full[(i-1)*DOWNSAMPLE+h, (j-1)*DOWNSAMPLE+h, (k-1)*DOWNSAMPLE+h]
end
phantom_full = nothing; GC.gc()
@info "  loaded + reversed + downsampled in $(round(time()-t0,digits=1))s  → $(size(phantom_labeled))"

# ── 4. Materials dict ─────────────────────────────────────────────────
# 4a. Existing 256..10355 cross-product bins from V1 voxelization (trees)
# 4b. NEW 10356..10456 chamber-only bins (100% blood, iodine ∈ [0,13])
const I_IODINE_EV = 491.0
const ZA_IODINE = 53.0 / 126.90447
const _MAT_BLOOD = BS.XA.Materials.blood
const _MAT_MUSCLE = BS.XA.Materials.muscle
function blend_bmi(f_blood::Real, C_iodine_mg_per_mL::Real)
    f = Float64(f_blood); C = Float64(C_iodine_mg_per_mL)
    m_b = f * BS.XA.val(_MAT_BLOOD.density)
    m_m = (1-f) * BS.XA.val(_MAT_MUSCLE.density)
    m_i = f * C * 1e-3
    total = m_b + m_m + m_i
    w_b = m_b/total; w_m = m_m/total; w_i = m_i/total
    ρ = total * Unitful.unit(_MAT_BLOOD.density)
    ZA = w_b * _MAT_BLOOD.ZA_ratio + w_m * _MAT_MUSCLE.ZA_ratio + w_i * ZA_IODINE
    Imix = (w_b * BS.XA.val(_MAT_BLOOD.I) + w_m * BS.XA.val(_MAT_MUSCLE.I) +
            w_i * I_IODINE_EV) * Unitful.unit(_MAT_BLOOD.I)
    comp = Dict{Int, Float64}()
    for el in union(keys(_MAT_BLOOD.composition), keys(_MAT_MUSCLE.composition), [53])
        cb = get(_MAT_BLOOD.composition, el, 0.0)
        cm = get(_MAT_MUSCLE.composition, el, 0.0)
        ci = (el == 53) ? 1.0 : 0.0
        comp[el] = w_b * cb + w_m * cm + w_i * ci
    end
    BS.XA.Material("bmi_$(round(100*f;digits=1))pc_$(round(C;digits=3))mgmL", ZA, Imix, ρ, comp)
end

mx = cfg["mixture_materials"]
const N_BLOOD_OLD = Int(mx["n_blood_bins"])
const N_IODINE_OLD = Int(mx["n_iodine_bins"])
const LABEL_BASE_OLD = Int(mx["label_base"])
const IODINE_MAX_OLD = Float64(mx["iodine_max_mg_per_mL"])

# new chamber bins
const CH_LABEL_BASE = 10356
const CH_N_BINS = 101                  # bin_i ∈ 0..100
const CH_IODINE_MAX = 13.0             # covers V2 RV peak ~12 mgI/mL

# materials dict
materials = Dict{Int, Any}()
for (k, v) in cfg["materials"]; materials[parse(Int, k)] = Symbol(v); end
for bin_b in 1:N_BLOOD_OLD, bin_i in 0:N_IODINE_OLD
    f = bin_b / Float64(N_BLOOD_OLD)
    C = bin_i / Float64(N_IODINE_OLD) * IODINE_MAX_OLD
    label = LABEL_BASE_OLD + (bin_b - 1) * (N_IODINE_OLD + 1) + bin_i
    materials[label] = blend_bmi(f, C)
end
for bin_i in 0:(CH_N_BINS - 1)
    C = bin_i / Float64(CH_N_BINS - 1) * CH_IODINE_MAX
    label = CH_LABEL_BASE + bin_i
    materials[label] = blend_bmi(1.0, C)  # 100% blood + iodine
end
v2_chamber_max = max(maximum(bae.C_aorta), maximum(bae.C_RV), maximum(bae.C_LV),
                     maximum(bae.C_pulm_artery), maximum(bae.C_pulm_vein))
@info "Materials" n_total=length(materials) old_bin_range="$(LABEL_BASE_OLD)..$(LABEL_BASE_OLD + N_BLOOD_OLD*(N_IODINE_OLD+1) - 1)" chamber_bin_range="$(CH_LABEL_BASE)..$(CH_LABEL_BASE+CH_N_BINS-1)" v2_chamber_max
@info "  (chamber iodine_max=$(CH_IODINE_MAX) covers V2 peak $(round(v2_chamber_max,digits=2)) ✓)"

# ── 5. Locate chamber voxels in downsampled phantom ───────────────────
const XCAT_TO_BAE = Dict{Symbol, Symbol}(
    :LV_blood_pool  => :LV,
    :RV_blood_pool  => :RV,
    :LA_blood_pool  => :pulm_vein,
    :RA_blood_pool  => :RV,
    :pulm_artery    => :pulm_artery,
    :pulm_veins     => :pulm_vein,
    :great_vessels  => :aorta,
)
chamber_records = cfg["chamber_patch"]
chamber_voxels = Dict{Symbol, Vector{Int}}()
inv_label_to_chamber = Dict{UInt16, Vector{Symbol}}()
# Multiple chambers can share the same init label (e.g. RV+RA both → 10342)
for (key, rec) in chamber_records
    startswith(key, "label_") || continue
    sym = Symbol(rec["name"])
    init_lbl = UInt16(rec["new_label"])
    chamber_voxels[sym] = Int[]
    if !haskey(inv_label_to_chamber, init_lbl)
        inv_label_to_chamber[init_lbl] = Symbol[]
    end
    push!(inv_label_to_chamber[init_lbl], sym)
end

@info "Indexing chamber voxels"
t0 = time()
# BUG to avoid: RA and RV share label 10342, so a single linear-index list per label
# would double-assign. Instead use the SOURCE XCAT label to disambiguate. Load XCAT
# UInt8 phantom and use that as the per-voxel chamber-ID source.
const XCAT_SRC_PATH = "/home/molloi-lab/smb_mount/shared_drive/Shu Nie/PVAT_Analysis/digital phantoms/vmale50_1600x1400x500_8bit_little_endian_act_1.raw"
const XCAT_SRC_DIMS = (1600, 1400, 500)
const XCAT_TO_CHAMBER_SYM = Dict{UInt8,Symbol}(
    UInt8(19) => :LV_blood_pool, UInt8(20) => :RV_blood_pool,
    UInt8(21) => :LA_blood_pool, UInt8(22) => :RA_blood_pool,
    UInt8(24) => :pulm_artery,   UInt8(25) => :pulm_veins,
    UInt8(28) => :great_vessels,
)
xcat_src = Array{UInt8}(undef, XCAT_SRC_DIMS...)
open(XCAT_SRC_PATH, "r") do io; read!(io, xcat_src); end
# Downsample XCAT same way to align with phantom_labeled coords (which were
# downsampled AFTER reverse(...; dims=(2,3))). For ROI mask building below we
# also downsample with the SAME reverse. Apply reverse(xcat_src; dims=(2,3)) too
# — but inspecting add_chambers_to_phantom shows the chamber_patched .raw is
# in *unreversed* XCAT coords (no reverse during patching), so we need to:
#   - Read XCAT source (no reverse), downsample.
#   - For phantom_labeled coords, we reversed (2,3). So XCAT must match THAT.
xcat_src = reverse(xcat_src; dims=(2,3))
xcat_ds = Array{UInt8}(undef, nx_ds, ny_ds, nz_ds)
@inbounds for k in 1:nz_ds, j in 1:ny_ds, i in 1:nx_ds
    xcat_ds[i,j,k] = xcat_src[(i-1)*DOWNSAMPLE+h, (j-1)*DOWNSAMPLE+h, (k-1)*DOWNSAMPLE+h]
end
xcat_src = nothing; GC.gc()

@inbounds for li in eachindex(phantom_labeled)
    v = xcat_ds[li]
    haskey(XCAT_TO_CHAMBER_SYM, v) || continue
    sym = XCAT_TO_CHAMBER_SYM[v]
    push!(chamber_voxels[sym], li)
end
for (sym, idxs) in chamber_voxels
    @info "  $(sym): $(length(idxs)) voxels"
end
@info "  indexing done in $(round(time()-t0,digits=1))s"

# ── 6. Per-time chamber relabel ───────────────────────────────────────
@inline function bin_i_for_C(C::Real)
    bi = round(Int, clamp(C / CH_IODINE_MAX, 0.0, 1.0) * (CH_N_BINS - 1))
    return clamp(bi, 0, CH_N_BINS - 1)
end
@inline function chamber_label_for_C(C::Real)
    UInt16(CH_LABEL_BASE + bin_i_for_C(C))
end

function patch_chambers_at_t!(phantom, bae, t)
    C_by_chamber = Dict{Symbol, Float64}()
    for (xcat_sym, bae_sym) in XCAT_TO_BAE
        C_field = bae_sym == :aorta       ? bae.C_aorta :
                  bae_sym == :RV          ? bae.C_RV :
                  bae_sym == :LV          ? bae.C_LV :
                  bae_sym == :pulm_artery ? bae.C_pulm_artery :
                  bae_sym == :pulm_vein   ? bae.C_pulm_vein :
                  error("unmapped Bae $bae_sym")
        C_by_chamber[xcat_sym] = c_at(bae.times, C_field, t)
    end
    for (sym, C) in C_by_chamber
        new_label = chamber_label_for_C(C)
        idxs = chamber_voxels[sym]
        @inbounds for li in idxs
            phantom[li] = new_label
        end
    end
    return C_by_chamber
end

# ── 7. ROI masks (in recon coords, from XCAT chamber labels) ──────────
@info "Building chamber ROI masks in recon coords"
const RECON_MATRIX = (512, 512, 250)
const RECON_FOV_CM = 35.0
const RECON_Z_CM   = 5.0

function build_roi_masks(recon_dims)
    # Use the downsampled xcat_ds we already built above (same coords as phantom_labeled)
    return _build_masks_inner(recon_dims, xcat_ds)
end
function _build_masks_inner(recon_dims, xcat_grid)
    nxg, nyg, nzg = size(xcat_grid)
    # CRITICAL: BS.Phantom interprets the downsampled (800,700,250) array as if
    # at full-res voxel pitch (0.02 cm), so its physical extent is only
    # (16, 14, 5) cm — half of the original 32×28×10 cm phantom. The other
    # half is "hidden" but the heart sits centered, so it remains in the visible
    # region. The MASK must use this same 0.02-cm interpretation, otherwise it
    # samples positions BS rendered as air → giving lung/air HU back instead of
    # chamber HU.
    bs_voxel = 0.02
    rvx = (RECON_FOV_CM/recon_dims[1], RECON_FOV_CM/recon_dims[2], RECON_Z_CM/recon_dims[3])
    org = (-nxg*bs_voxel/2, -nyg*bs_voxel/2, -nzg*bs_voxel/2)
    masks = Dict{Symbol, BitArray{3}}()
    counts = Dict{Symbol, Int}()
    for n in values(XCAT_TO_CHAMBER_SYM); masks[n] = falses(recon_dims...); counts[n] = 0; end
    @inbounds for kr in 1:recon_dims[3], jr in 1:recon_dims[2], ir in 1:recon_dims[1]
        x = -RECON_FOV_CM/2 + (ir-0.5)*rvx[1]
        y = -RECON_FOV_CM/2 + (jr-0.5)*rvx[2]
        z = -RECON_Z_CM/2  + (kr-0.5)*rvx[3]
        i_ds = round(Int, (x - org[1])/bs_voxel + 0.5)
        j_ds = round(Int, (y - org[2])/bs_voxel + 0.5)
        k_ds = round(Int, (z - org[3])/bs_voxel + 0.5)
        (1<=i_ds<=nxg && 1<=j_ds<=nyg && 1<=k_ds<=nzg) || continue
        v = xcat_grid[i_ds, j_ds, k_ds]
        haskey(XCAT_TO_CHAMBER_SYM, v) || continue
        name = XCAT_TO_CHAMBER_SYM[v]
        masks[name][ir,jr,kr] = true
        counts[name] += 1
    end
    @info "ROI masks built" counts...
    return masks
end
roi_masks = build_roi_masks(RECON_MATRIX)
# DA-only ROI: take great_vessels label restricted to lower z (below heart).
# Heart sits in mid-z; XCAT label 28 covers SVC + DA + IVC region. We isolate DA
# by taking the BOTTOM 30% of label-28 voxels in z (closer to diaphragm) as the
# DA proxy. This is a common practical workaround when label 28 lumps vessels.
let
    gv = roi_masks[:great_vessels]
    z_positions = Int[]
    nz = size(gv, 3)
    for k in 1:nz; any(gv[:,:,k]) && push!(z_positions, k); end
    if !isempty(z_positions)
        z_lo, z_hi = z_positions[1], z_positions[end]
        z_cut = z_lo + Int(round((z_hi - z_lo) * 0.30))
        da_mask = falses(size(gv)...)
        @inbounds for k in z_lo:z_cut, j in 1:size(gv,2), i in 1:size(gv,1)
            da_mask[i,j,k] = gv[i,j,k]
        end
        roi_masks[:DA] = da_mask
        @info "  DA mask (z $(z_lo)..$(z_cut)): $(count(da_mask)) voxels"
    end
end

# ── 8. GPU + scanner ───────────────────────────────────────────────────
const GPU_BACKEND = let
    detected = (name="CPU", to_gpu=identity)
    for (pkg, uuid, ctor) in (
        (:CUDA, "052768ef-5323-5732-b1bb-66c8b64840ba", :CuArray),
        (:Metal, "dde4c033-4e86-420c-a63e-0dd931031962", :MtlArray),
    )
        pkg_id = Base.PkgId(Base.UUID(uuid), String(pkg))
        Base.locate_package(pkg_id) === nothing && continue
        try
            m = Base.require(pkg_id)
            if Base.invokelatest(getfield(m, :functional))
                detected = (name=string(pkg), to_gpu=getfield(m, ctor))
                break
            end
        catch
        end
    end
    detected
end
to_gpu(x) = GPU_BACKEND.to_gpu(x)
@info "GPU backend: $(GPU_BACKEND.name)"

scanner = BS.Scanner(
    source_to_isocenter = 625.6, source_to_detector = 1100.0,
    detector_rows = 256, detector_cols = 834,
    detector_row_size = 0.625, detector_col_size = 0.6,
    detector_shape = BS.CURVED_DETECTOR,
    focal_spot_width = 1.0, focal_spot_length = 1.0, target_angle = 10.0,
    flat_filter_material = :aluminum, flat_filter_thickness = 2.5,
    bowtie_filter = :ge_revolution_large,
    detector_material = :lumex, detector_depth = 3.0,
    fill_factor_row = 0.9, fill_factor_col = 0.9,
    electronic_noise = 0, detection_gain = 10.0,
)
ct_protocol = BS.CTProtocol(kVp = KVP, mA = MA, views = VIEWS, rotation_time = 1.0,
                             collimation_mm = 50.0, additional_filters = [("Al", 4.5)])
sim_opts = BS.SimOptions(fidelity = :eict, seed = 1234)
recon_opts = BS.ReconOptions(matrix_size = RECON_MATRIX, fov_cm = RECON_FOV_CM, z_cm = RECON_Z_CM)

# ── 9. Sweep ──────────────────────────────────────────────────────────
csv_path = joinpath(A.output_dir, "chamber_hu_timeseries.csv")
chambers_ordered = [:DA, :great_vessels, :LV_blood_pool, :RV_blood_pool,
                     :LA_blood_pool, :RA_blood_pool, :pulm_artery, :pulm_veins]
open(csv_path, "w") do io
    header = "t_s," * join(["HU_$(c)" for c in chambers_ordered], ",") *
             "," * join(["bae_C_$(c)_mgI_ml" for c in chambers_ordered], ",")
    println(io, header)
end

@info "Starting time sweep" n_iters=length(SCAN_TIMES) dt=A.dt t_end=A.t_end
bhc_model = nothing
total_t0 = time()

for (idx, t) in enumerate(SCAN_TIMES)
    global bhc_model  # Julia 1.11+ soft-scope: explicitly bind to module global
    iter_t0 = time()
    @info "iter $(idx)/$(length(SCAN_TIMES)) — t = $(round(t, digits=2)) s"

    # 9a. patch chambers
    C_by_chamber = patch_chambers_at_t!(phantom_labeled, bae, t)
    @info "  chambers: " * join(["$(k)=$(round(v;digits=2))" for (k,v) in C_by_chamber], " ")

    # 9b. forward project
    phantom_bs = BS.Phantom(to_gpu(phantom_labeled), materials, VOXEL_SIZE_CM)
    ws = BS.create_eict_workspace(scanner, ct_protocol, sim_opts, recon_opts, phantom_bs)
    BS.simulate!(ws, phantom_bs, ct_protocol, sim_opts)
    sino_cpu = Array(ws.sinogram); geom = ws.geom
    ws = nothing; phantom_bs = nothing; GC.gc(true)

    # 9c. BHC calibration (first iter only)
    if bhc_model === nothing
        prot_for_bhc = BS.CTProtocol(kVp = KVP, additional_filters = [("Al", 4.5)])
        bhc_model = BS.calibrate_bhc_two_material(
            sim_opts, prot_for_bhc; scanner=scanner, geom=geom,
            order=2, hu_low=450.0, hu_high=600.0)
        @info "  BHC calibrated  ref_E_keV=$(round(bhc_model.reference_energy_keV;digits=1))"
    end

    # 9d. BHC + FDK
    sino_gpu = to_gpu(sino_cpu)
    sino_bhc = BS.apply_bhc_two_material(sino_gpu, bhc_model, geom, RECON_MATRIX)
    sino_gpu = to_gpu(sino_bhc)
    ws_fdk = BS.create_fdk_recon_workspace(sino_gpu, geom, RECON_MATRIX)
    recon_μ = BS.reconstruct!(ws_fdk, sino_gpu, geom)
    BS.apply_bhc_image_domain(recon_μ, geom, RECON_MATRIX, bhc_model.μ_water_ref;
                               hu_low=50.0, hu_high=150.0, scale_factor=0.2)
    hu = Float32.(BS.to_hounsfield(Array(recon_μ); μ_water = bhc_model.μ_water_ref))
    BS.add_system_noise_floor!(hu, 28.0; seed=1234)
    BS.apply_radial_cupping_correction!(hu; fov_cm=RECON_FOV_CM)
    ws_fdk = nothing; sino_gpu = nothing; recon_μ = nothing
    sino_cpu = nothing; sino_bhc = nothing
    GC.gc(true)

    # 9e. sample ROIs
    row_hu = Float64[]
    for ch in chambers_ordered
        if !haskey(roi_masks, ch)
            push!(row_hu, NaN); continue
        end
        mask = roi_masks[ch]
        v = any(mask) ? mean(hu[mask]) : NaN
        push!(row_hu, v)
    end
    bae_C_row = Float64[]
    for ch in chambers_ordered
        if ch == :DA || ch == :great_vessels
            push!(bae_C_row, c_at(bae.times, bae.C_aorta, t))
        elseif ch == :LV_blood_pool
            push!(bae_C_row, c_at(bae.times, bae.C_LV, t))
        elseif ch == :RV_blood_pool || ch == :RA_blood_pool
            push!(bae_C_row, c_at(bae.times, bae.C_RV, t))
        elseif ch == :LA_blood_pool || ch == :pulm_veins
            push!(bae_C_row, c_at(bae.times, bae.C_pulm_vein, t))
        elseif ch == :pulm_artery
            push!(bae_C_row, c_at(bae.times, bae.C_pulm_artery, t))
        else
            push!(bae_C_row, NaN)
        end
    end
    open(csv_path, "a") do io
        cells = [@sprintf("%.3f", t)]
        for v in row_hu;   push!(cells, @sprintf("%.2f", v)); end
        for v in bae_C_row; push!(cells, @sprintf("%.4f", v)); end
        println(io, join(cells, ","))
    end
    @info "  HU @ t=$(round(t,digits=1))s: DA=$(round(Int,row_hu[1])) GV=$(round(Int,row_hu[2])) LV=$(round(Int,row_hu[3])) RV=$(round(Int,row_hu[4]))  ($(round(time()-iter_t0;digits=1))s)"

    # 9f. final-iter DICOM
    if idx == length(SCAN_TIMES)
        @info "  Final iter — writing DICOM"
        voxel_mm = (RECON_FOV_CM * 10 / RECON_MATRIX[2],
                    RECON_FOV_CM * 10 / RECON_MATRIX[1],
                    RECON_Z_CM * 10 / RECON_MATRIX[3])
        write_hu_to_dicom_series(
            joinpath(A.output_dir, "dicom_final_t$(round(Int, t))s"), hu, voxel_mm;
            study_uid=gen_uid(), series_uid=gen_uid(), frame_of_ref_uid=gen_uid(),
            series_description="UCI V2 time-series final t=$(round(Int,t))s",
            series_number=1, kvp=KVP, mA=MA,
            window_center=250.0, window_width=500.0, acquisition_dt=now(),
        )
        @info "  DICOM written"
    end
    hu = nothing; GC.gc(true)
end
@info "Sweep done in $(round((time()-total_t0)/60;digits=1)) min  CSV: $(csv_path)"
println("Done.")
