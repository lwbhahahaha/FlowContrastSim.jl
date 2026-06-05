#!/usr/bin/env julia
# analyze_vmi_outputs.jl — sample chamber HU from VMI reconstructions
# at multiple keV and the single-kVp 120 kVp baseline, then write a
# table comparing Bae predictions, 120 kVp polychromatic recon, and
# 70 keV VMI (the "monoenergetic equivalent" of 120 kVp).
#
# Diagnostic:
#   - If LV HU @ 70 keV VMI ≈ Bae prediction but 120 kVp polychromatic
#     is lower → BHC was the issue (single-kVp BHC residuals).
#   - If LV HU @ 70 keV VMI ≈ 120 kVp polychromatic → BHC is OK, the
#     gap vs Bae is upstream (chamber-patching encoding or flow physics).
#
# CLI:
#   julia analyze_vmi_outputs.jl   RATE   [cta_out_dual_dir]
#
# Default cta_out_dual_dir = phantom_ct_input/bae_uci_r{rate}/cta_out_dual

using TOML
using Statistics
using Printf

const ROOT = abspath(joinpath(@__DIR__, "..", ".."))

length(ARGS) >= 1 || error("Usage: julia analyze_vmi_outputs.jl  RATE")
const RATE = ARGS[1]
const CTA_DUAL_DIR = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(ROOT, "phantom_ct_input", "bae_uci_r$(RATE)", "cta_out_dual")
const CTA_120_DIR = joinpath(ROOT, "phantom_ct_input", "bae_uci_r$(RATE)", "cta_out")
const PEAK_DIR = joinpath(ROOT, "phantom_ct_input", "bae_uci_r$(RATE)", "peak_iodine")

const XCAT_CHAMBER_LABELS = Dict{UInt8,Symbol}(
    UInt8(19) => :LV_blood_pool,
    UInt8(20) => :RV_blood_pool,
    UInt8(21) => :LA_blood_pool,
    UInt8(22) => :RA_blood_pool,
    UInt8(24) => :pulm_artery,
    UInt8(25) => :pulm_veins,
    UInt8(28) => :great_vessels,
)

const XCAT_SRC_PATH = "/home/molloi-lab/smb_mount/shared_drive/Shu Nie/PVAT_Analysis/digital phantoms/vmale50_1600x1400x500_8bit_little_endian_act_1.raw"
const XCAT_SRC_DIMS = (1600, 1400, 500)

function load_hu(path::String, dims::NTuple{3,Int})
    arr = Array{Float32}(undef, dims...)
    open(path, "r") do io; read!(io, arr); end
    arr
end

# Chamber mask in recon coords (same logic as analyze_uci_outputs.jl)
function chamber_masks(recon_dims, bs_voxel_cm=0.02, recon_fov_cm=35.0, recon_z_cm=5.0)
    @info "Loading XCAT source + building masks…"
    src = Array{UInt8}(undef, XCAT_SRC_DIMS...)
    open(XCAT_SRC_PATH, "r") do io; read!(io, src); end
    nx_src, ny_src, nz_src = XCAT_SRC_DIMS
    ds = 2
    nx_ds = nx_src ÷ ds; ny_ds = ny_src ÷ ds; nz_ds = nz_src ÷ ds
    h = ds ÷ 2 + 1
    grid = Array{UInt8}(undef, nx_ds, ny_ds, nz_ds)
    @inbounds for k in 1:nz_ds, j in 1:ny_ds, i in 1:nx_ds
        i_s = (i-1)*ds + h
        j_s = ny_src - ((j-1)*ds + h) + 1
        k_s = nz_src - ((k-1)*ds + h) + 1
        grid[i,j,k] = src[i_s, j_s, k_s]
    end
    src = nothing; GC.gc()
    nrx, nry, nrz = recon_dims
    rvx = (recon_fov_cm/nrx, recon_fov_cm/nry, recon_z_cm/nrz)
    ext = (nx_ds*bs_voxel_cm, ny_ds*bs_voxel_cm, nz_ds*bs_voxel_cm)
    org = (-ext[1]/2, -ext[2]/2, -ext[3]/2)
    masks = Dict{Symbol, BitArray{3}}()
    counts = Dict{Symbol, Int}()
    for n in values(XCAT_CHAMBER_LABELS); masks[n] = falses(nrx, nry, nrz); counts[n] = 0; end
    @inbounds for kr in 1:nrz, jr in 1:nry, ir in 1:nrx
        x = -recon_fov_cm/2 + (ir-0.5)*rvx[1]
        y = -recon_fov_cm/2 + (jr-0.5)*rvx[2]
        z = -recon_z_cm/2 + (kr-0.5)*rvx[3]
        i_ds = round(Int, (x - org[1]) / bs_voxel_cm + 0.5)
        j_ds = round(Int, (y - org[2]) / bs_voxel_cm + 0.5)
        k_ds = round(Int, (z - org[3]) / bs_voxel_cm + 0.5)
        (1<=i_ds<=nx_ds && 1<=j_ds<=ny_ds && 1<=k_ds<=nz_ds) || continue
        v = grid[i_ds, j_ds, k_ds]
        if haskey(XCAT_CHAMBER_LABELS, v)
            name = XCAT_CHAMBER_LABELS[v]
            masks[name][ir, jr, kr] = true
            counts[name] += 1
        end
    end
    @info "Mask voxel counts" counts...
    return masks
end

function roi_stats(hu, masks)
    out = Dict{Symbol, NamedTuple}()
    for (name, mask) in masks
        if any(mask)
            v = hu[mask]
            out[name] = (n=length(v), mean=mean(v), std=std(v),
                         min=minimum(v), max=maximum(v))
        else
            out[name] = (n=0, mean=NaN, std=NaN, min=NaN, max=NaN)
        end
    end
    return out
end

# ── Load all recons ────────────────────────────────────────────────
@info "Reading recon_meta from $(CTA_DUAL_DIR)"
meta_dual = TOML.parsefile(joinpath(CTA_DUAL_DIR, "recon_meta.toml"))
recon_dims = Tuple(Int.(meta_dual["recon"]["shape"]))
keVs = Float64.(meta_dual["vmi"]["keVs"])

vmi = Dict{Float64, Array{Float32,3}}()
for E in keVs
    path = joinpath(CTA_DUAL_DIR, "vmi_HU_$(Int(E))keV_f32.raw")
    vmi[E] = load_hu(path, recon_dims)
    @info "  loaded VMI $(Int(E)) keV  $(round(filesize(path)/1024^2;digits=0)) MB"
end

# Single-kVp 120 baseline
hu_fbp120 = load_hu(joinpath(CTA_120_DIR, "recon_fbp_hu_f32.raw"), recon_dims)
hu_hir120 = load_hu(joinpath(CTA_120_DIR, "recon_hir_hu_f32.raw"), recon_dims)
@info "loaded 120 kVp poly FBP + HIR"

# ── Build masks ────────────────────────────────────────────────────
masks = chamber_masks(recon_dims)

# ── Bae chamber predictions ─────────────────────────────────────────
pm = TOML.parsefile(joinpath(PEAK_DIR, "peak_metadata.toml"))
cc = pm["chamber_concentrations"]
# Predicted HU at each VMI energy:
# At 70 keV (≈ 120 kVp polychromatic effective E), slope ≈ 25 HU/(mgI/mL)
# Other keV need their own slope from the iodine + water μ tables.
# Approximation: HU(E) = 1000 * C * (μ/ρ)_I(E) / (μ/ρ)_water(E)  (low-conc limit, no baseline blood)
# Actually I'll just compute the slope from typical literature values.
# 40 keV: ~52 ;  55 keV: ~36 ; 70 keV: ~26 ; 85 keV: ~22 ; 100 keV: ~19 ; 140 keV: ~15
const HU_PER_MGI_AT_KEV = Dict(
    40.0 => 52.0, 55.0 => 36.0, 70.0 => 26.0, 85.0 => 22.0,
    100.0 => 19.0, 140.0 => 15.0,
)

# ── Sample everything ──────────────────────────────────────────────
@info "Sampling ROIs…"
stats_120fbp = roi_stats(hu_fbp120, masks)
stats_120hir = roi_stats(hu_hir120, masks)
stats_vmi = Dict(E => roi_stats(vmi[E], masks) for E in keVs)

# ── Output ─────────────────────────────────────────────────────────
out_md = joinpath(ROOT, "phantom_ct_input", "vmi_analysis_r$(RATE).md")
open(out_md, "w") do io
    println(io, "# VMI diagnostic — rate=$(RATE) mL/s, 100 kg patient, UCI triphasic\n")
    println(io, "Compares: Bae prediction · 120 kVp poly (FBP + HIR) · DECT VMI at multiple keV")
    println(io)
    println(io, "## Per-chamber HU (mean over recon-coord mask)\n")
    chambers = [:LV_blood_pool, :RV_blood_pool, :LA_blood_pool, :RA_blood_pool,
                :pulm_artery, :pulm_veins, :great_vessels]
    bae_key = Dict(:LV_blood_pool=>"left_heart_mgI_ml", :RV_blood_pool=>"right_heart_mgI_ml",
                   :LA_blood_pool=>"pulm_vein_mgI_ml", :RA_blood_pool=>"right_heart_mgI_ml",
                   :pulm_artery=>"pulm_artery_mgI_ml", :pulm_veins=>"pulm_vein_mgI_ml",
                   :great_vessels=>"aorta_root_mgI_ml")

    # Header row
    cols = ["Chamber", "Bae C (mgI/mL)", "Bae HU₀ (120)", "120 FBP", "120 HIR"]
    for E in keVs; push!(cols, "VMI $(Int(E)) keV"); end
    println(io, "| " * join(cols, " | ") * " |")
    println(io, "|" * repeat(" --- |", length(cols)))

    for ch in chambers
        c = Float64(cc[bae_key[ch]])
        bae_hu120 = 40 + 25 * c
        row = [String(ch), @sprintf("%.2f", c), @sprintf("%.0f", bae_hu120),
               @sprintf("%.0f", stats_120fbp[ch].mean),
               @sprintf("%.0f", stats_120hir[ch].mean)]
        for E in keVs
            push!(row, @sprintf("%.0f", stats_vmi[E][ch].mean))
        end
        println(io, "| " * join(row, " | ") * " |")
    end

    # Predicted VMI HU per chamber per keV (for reference)
    println(io)
    println(io, "## Bae-predicted VMI HU per chamber per keV (literature slopes)\n")
    cols2 = ["Chamber", "C (mgI/mL)"]
    for E in keVs; push!(cols2, "$(Int(E)) keV"); end
    println(io, "| " * join(cols2, " | ") * " |")
    println(io, "|" * repeat(" --- |", length(cols2)))
    for ch in chambers
        c = Float64(cc[bae_key[ch]])
        row = [String(ch), @sprintf("%.2f", c)]
        for E in keVs
            slope = get(HU_PER_MGI_AT_KEV, E, 25.0)
            pred = 40 + slope * c
            push!(row, @sprintf("%.0f", pred))
        end
        println(io, "| " * join(row, " | ") * " |")
    end
    println(io)
    println(io, "## Diagnostic\n")
    lv_bae120 = 40 + 25 * Float64(cc["left_heart_mgI_ml"])
    lv_120hir = stats_120hir[:LV_blood_pool].mean
    lv_70keV  = stats_vmi[70.0][:LV_blood_pool].mean
    @printf(io, "LV HU:  Bae pred @ 120 kVp ≈ %.0f  ·  120 kVp HIR (poly) = %.0f  ·  70 keV VMI = %.0f\n",
            lv_bae120, lv_120hir, lv_70keV)
    gap_120 = lv_bae120 - lv_120hir
    gap_70 = lv_bae120 - lv_70keV
    @printf(io, "  Δ(Bae − 120 HIR) = %.0f HU\n", gap_120)
    @printf(io, "  Δ(Bae − 70 keV VMI) = %.0f HU\n", gap_70)
    println(io)
    if abs(gap_70) < abs(gap_120) - 10
        println(io, "→ VMI 70 keV closer to Bae than 120 kVp poly. **BHC residual is the dominant cause** of the 120 kVp gap.")
    elseif abs(gap_70) > abs(gap_120) + 5
        println(io, "→ VMI 70 keV further from Bae. Suggests decomp / VMI bias rather than 120 kVp BHC.")
    else
        println(io, "→ VMI 70 keV ≈ 120 kVp poly HU. The 120 kVp gap is NOT BHC — it is upstream (chamber-patching / Bae).")
    end
end
print(read(out_md, String))
@info "Wrote $out_md"
