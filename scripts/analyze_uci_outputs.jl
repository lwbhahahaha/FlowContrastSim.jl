#!/usr/bin/env julia
# analyze_uci_outputs.jl — load both rate=2.5 and rate=1.0 reconstructed
# HU volumes, sample chamber ROIs, render side-by-side PNG and write a
# comparison table.
#
# Sampling strategy:
#   - Load the chamber-patched UInt16 phantom (source coord, 1600×1400×500)
#     and the reconstructed HU Float32 volume (512×512×250 after FOV=35cm
#     reconstruction at 0.68×0.68×0.20 mm voxels).
#   - For each chamber label (LV blood pool 19/cross-product → 10347 etc),
#     find the phantom-coordinate centroid + extent.
#   - The recon does:  reverse(phantom; dims=(2,3))  then downsample by 2.
#     So recon[i, j, k] ≈ source[(2i-1), (ny - 2j + 2), (nz - 2k + 2)] for
#     the simple case (with recon FOV cropping).
#   - Cleaner: just iterate recon voxels, compute back-mapped source coord,
#     find chamber label, accumulate HU into per-chamber bins.
#
# Output:
#   phantom_ct_input/uci_compare.md       — markdown table
#   phantom_ct_input/uci_compare.png      — FBP and HIR mid-slices, both rates

using TOML
using Printf
using Statistics

const ROOT = abspath(joinpath(@__DIR__, "..", ".."))

# Cross-product chamber label remap used by add_chambers_to_phantom.jl
# (these are the *new* UInt16 labels after patching; XCAT labels 19/20/21/22/24/25/28
# all get cross-product encoding with bin_b=100 and a per-chamber bin_i.)
# Read directly from the patched manifest to be robust to label drift.

function load_chamber_labels(manifest_path::String)
    cfg = TOML.parsefile(manifest_path)
    haskey(cfg, "chamber_patch") || error("manifest missing [chamber_patch] section: $manifest_path")
    cp = cfg["chamber_patch"]
    # Each entry: label_19 = { name = "LV_blood_pool", new_label = 10347, ... }
    labels = Dict{Symbol, Int}()
    for (k, v) in cp
        startswith(k, "label_") || continue
        name = Symbol(v["name"])
        labels[name] = Int(v["new_label"])
    end
    return labels
end

function load_phantom_u16(path::String, dims::NTuple{3,Int})
    nx, ny, nz = dims
    arr = Array{UInt16}(undef, nx, ny, nz)
    open(path, "r") do io; read!(io, arr); end
    return arr
end

function load_hu_f32(path::String, dims::NTuple{3,Int})
    nx, ny, nz = dims
    arr = Array{Float32}(undef, nx, ny, nz)
    open(path, "r") do io; read!(io, arr); end
    return arr
end

# Build a chamber mask in *recon coordinates*. Recon was built by:
#   1. read source phantom UInt16 (nx_src × ny_src × nz_src)
#   2. reverse along dims (2, 3)
#   3. downsample by factor=2  →  mask grid (nx_src÷2 × ny_src÷2 × nz_src÷2)
#   4. project onto 512×512×nz grid at FOV=35 cm
# Step 4 is just a coordinate remap; for the centre 14 cm of FOV the
# downsampled mask grid (200 µm voxels, 16 cm extent) sits inside the
# recon FOV (35 cm at 0.68 mm voxel = 35×35 cm). The mask is centred in
# the recon FOV by phantom.origin.
#
# Simplified: build chamber mask at the *downsampled-and-reversed* grid
# (800 × 700 × 250 = recon's z-matched mask), then map to recon by
# linear isocenter-based mapping.

function chamber_mask_in_recon_coords(phantom_src::Array{UInt8,3},
                                       xcat_chamber_labels::Dict{UInt8,Symbol},
                                       recon_dims::NTuple{3,Int};
                                       ds_factor::Int = 2,
                                       bs_voxel_cm::Float64 = 0.02,
                                       recon_fov_cm::Float64 = 35.0,
                                       recon_z_cm::Float64 = 5.0)
    # Uses the ORIGINAL XCAT phantom (UInt8, 1600×1400×500) which has
    # distinct labels 19,20,21,22,24,25,28 — robust against chamber labels
    # being aliased to the same cross-product code after add_chambers patch.
    nx_src, ny_src, nz_src = size(phantom_src)
    @info "Building chamber masks in recon coords from XCAT UInt8 source"
    # BS Phantom kept source voxel size (0.02 cm) even after downsample — its
    # physical extent is (nx_ds × 0.02, ny_ds × 0.02, nz_ds × 0.02) cm,
    # centered at isocenter. The recon FOV (35 cm) is much larger than this
    # extent; most recon voxels fall OUTSIDE the phantom (read as air).
    nx_ds = nx_src ÷ ds_factor
    ny_ds = ny_src ÷ ds_factor
    nz_ds = nz_src ÷ ds_factor
    h = ds_factor ÷ 2 + 1
    ds_grid = Array{UInt8}(undef, nx_ds, ny_ds, nz_ds)
    @inbounds for k in 1:nz_ds, j in 1:ny_ds, i in 1:nx_ds
        i_s = (i - 1) * ds_factor + h
        j_s = ny_src - ((j - 1) * ds_factor + h) + 1   # reverse(y)
        k_s = nz_src - ((k - 1) * ds_factor + h) + 1   # reverse(z)
        ds_grid[i, j, k] = phantom_src[i_s, j_s, k_s]
    end

    nrx, nry, nrz = recon_dims
    recon_voxel_cm = (recon_fov_cm / nrx, recon_fov_cm / nry, recon_z_cm / nrz)
    ds_extent_cm = (nx_ds * bs_voxel_cm, ny_ds * bs_voxel_cm, nz_ds * bs_voxel_cm)
    ds_origin_cm = (-ds_extent_cm[1]/2, -ds_extent_cm[2]/2, -ds_extent_cm[3]/2)

    masks = Dict{Symbol, BitArray{3}}()
    for name in values(xcat_chamber_labels)
        masks[name] = falses(nrx, nry, nrz)
    end
    n_inside = Dict{Symbol, Int}(name => 0 for name in values(xcat_chamber_labels))

    @inbounds for kr in 1:nrz, jr in 1:nry, ir in 1:nrx
        x_cm = -recon_fov_cm/2 + (ir - 0.5) * recon_voxel_cm[1]
        y_cm = -recon_fov_cm/2 + (jr - 0.5) * recon_voxel_cm[2]
        z_cm = -recon_z_cm/2 + (kr - 0.5) * recon_voxel_cm[3]
        i_ds = round(Int, (x_cm - ds_origin_cm[1]) / bs_voxel_cm + 0.5)
        j_ds = round(Int, (y_cm - ds_origin_cm[2]) / bs_voxel_cm + 0.5)
        k_ds = round(Int, (z_cm - ds_origin_cm[3]) / bs_voxel_cm + 0.5)
        (1 <= i_ds <= nx_ds && 1 <= j_ds <= ny_ds && 1 <= k_ds <= nz_ds) || continue
        v = ds_grid[i_ds, j_ds, k_ds]
        if haskey(xcat_chamber_labels, v)
            name = xcat_chamber_labels[v]
            masks[name][ir, jr, kr] = true
            n_inside[name] += 1
        end
    end
    @info "Chamber voxels in recon coords:" n_inside...
    return masks
end

function chamber_hu_stats(hu::Array{Float32,3}, masks::Dict{Symbol, BitArray{3}})
    out = Dict{Symbol, NamedTuple}()
    for (name, mask) in masks
        if any(mask)
            v = hu[mask]
            out[name] = (n = length(v), mean = mean(v), median = median(v),
                          std = std(v), min = minimum(v), max = maximum(v))
        else
            out[name] = (n = 0, mean = NaN, median = NaN, std = NaN, min = NaN, max = NaN)
        end
    end
    return out
end

const XCAT_CHAMBER_LABELS = Dict{UInt8,Symbol}(
    UInt8(19) => :LV_blood_pool,
    UInt8(20) => :RV_blood_pool,
    UInt8(21) => :LA_blood_pool,
    UInt8(22) => :RA_blood_pool,
    UInt8(24) => :pulm_artery,
    UInt8(25) => :pulm_veins,
    UInt8(28) => :great_vessels,
)

# Load original 1600×1400×500 UInt8 XCAT (chambers still distinct labels).
const XCAT_SRC_PATH = "/home/molloi-lab/smb_mount/shared_drive/Shu Nie/PVAT_Analysis/digital phantoms/vmale50_1600x1400x500_8bit_little_endian_act_1.raw"
const XCAT_SRC_DIMS = (1600, 1400, 500)
const XCAT_SRC = let
    arr = Array{UInt8}(undef, XCAT_SRC_DIMS...)
    open(XCAT_SRC_PATH, "r") do io; read!(io, arr); end
    arr
end
@info "Loaded XCAT UInt8 source: $(size(XCAT_SRC))"

function analyze_run(rate::String)
    work = joinpath(ROOT, "phantom_ct_input", "bae_uci_r$(rate)")
    cta_dir = joinpath(work, "cta_out")
    pm_dir = joinpath(work, "peak_iodine")

    recon_meta = TOML.parsefile(joinpath(cta_dir, "recon_meta.toml"))
    recon_dims = Tuple(Int.(recon_meta["recon"]["shape"]))
    hu_fbp = load_hu_f32(joinpath(cta_dir, "recon_fbp_hu_f32.raw"), recon_dims)
    hu_hir = load_hu_f32(joinpath(cta_dir, "recon_hir_hu_f32.raw"), recon_dims)

    @info "rate=$(rate): computing chamber masks in recon coords"
    masks = chamber_mask_in_recon_coords(XCAT_SRC, XCAT_CHAMBER_LABELS, recon_dims;
                                          bs_voxel_cm=0.02,
                                          recon_z_cm=Float64(recon_meta["recon"]["z_cm"]))

    # Bae chamber concentrations from peak_metadata
    pm = TOML.parsefile(joinpath(pm_dir, "peak_metadata.toml"))
    cc = pm["chamber_concentrations"]

    return (rate=rate, recon_meta=recon_meta, cc=cc,
            stats_fbp = chamber_hu_stats(hu_fbp, masks),
            stats_hir = chamber_hu_stats(hu_hir, masks),
            hu_fbp = hu_fbp, hu_hir = hu_hir)
end

rates = ["1.5", "2.0"]
@info "Analyzing rates: $rates"
results = [analyze_run(r) for r in rates]

# ── Markdown table ───────────────────────────────────────────────────
out_md = joinpath(ROOT, "phantom_ct_input", "uci_compare.md")
open(out_md, "w") do io
    println(io, "# UCI Cardiac CTP — phase-2 rate comparison")
    println(io)
    println(io, "Bae 1998 PBPK + Taylor-Aris PDE + V4 cross-product voxelization + BasisSim 120 kVp")
    println(io)
    println(io, "## Bae chamber predictions at scan time")
    println(io)
    r0, r1 = rates[1], rates[2]
    println(io, "| chamber | C @ rate=$(r0) (mgI/mL) | HU₀ pred | C @ rate=$(r1) (mgI/mL) | HU₀ pred |")
    println(io, "| --- | :---: | :---: | :---: | :---: |")
    for k in ("aorta_root_mgI_ml", "left_heart_mgI_ml", "right_heart_mgI_ml",
              "pulm_artery_mgI_ml", "pulm_vein_mgI_ml")
        nm = replace(replace(k, "_mgI_ml"=>""), "_"=>" ")
        c0 = Float64(results[1].cc[k]); h0 = 40 + 25 * c0
        c1 = Float64(results[2].cc[k]); h1 = 40 + 25 * c1
        @printf(io, "| %-14s | %.2f | %.0f | %.2f | %.0f |\n", nm, c0, h0, c1, h1)
    end
    println(io)
    for (i, recon) in ((1, "FBP"), (2, "HIR"))
        stats_field = recon == "FBP" ? :stats_fbp : :stats_hir
        println(io, "## Reconstructed HU stats — $(recon)")
        println(io)
        println(io, "| chamber | rate=$(r0) mean | rate=$(r0) std | rate=$(r0) n | rate=$(r1) mean | rate=$(r1) std | rate=$(r1) n |")
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
        for k in [:LV_blood_pool, :RV_blood_pool, :LA_blood_pool, :RA_blood_pool,
                   :pulm_artery, :pulm_veins, :great_vessels]
            s0 = getproperty(results[1], stats_field)[k]
            s1 = getproperty(results[2], stats_field)[k]
            @printf(io, "| %-15s | %6.1f | %5.1f | %6d | %6.1f | %5.1f | %6d |\n",
                    String(k), s0.mean, s0.std, s0.n, s1.mean, s1.std, s1.n)
        end
        println(io)
    end
    println(io, "## Beam-hardening criterion: HU_RV > 300")
    println(io)
    for (r, recon) in ((r0,"FBP"), (r0,"HIR"), (r1,"FBP"), (r1,"HIR"))
        idx = r == r0 ? 1 : 2
        stats = recon == "FBP" ? results[idx].stats_fbp : results[idx].stats_hir
        rv = stats[:RV_blood_pool]
        bh = rv.mean > 300 ? "⚠ BH RISK" : "✓ OK"
        @printf(io, "  rate=%s  %s:  RV mean = %.0f HU   %s\n", r, recon, rv.mean, bh)
    end
end
@info "Wrote $out_md"
print(read(out_md, String))
