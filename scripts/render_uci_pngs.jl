#!/usr/bin/env julia
# render_uci_pngs.jl — render mid-slice side-by-side PNG comparison for
# the two UCI phase-2 rate runs. Outputs phantom_ct_input/uci_compare.png.
#
# Uses CairoMakie (already a dep of phantom_ct_input/run_basis_sim).

import Pkg
const ROOT = abspath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "phantom_ct_input", "run_basis_sim"))

import CairoMakie as CM
using TOML
using Statistics
using Printf

const RATES = ("1.5", "2.0")

function load_recon(rate)
    work = joinpath(ROOT, "phantom_ct_input", "bae_uci_r$(rate)", "cta_out")
    rm = TOML.parsefile(joinpath(work, "recon_meta.toml"))
    shape = Tuple(Int.(rm["recon"]["shape"]))
    n = prod(shape)
    function load(name)
        data = Vector{Float32}(undef, n)
        open(joinpath(work, name), "r") do io; read!(io, data); end
        reshape(data, shape)
    end
    return (rate=rate, fbp=load("recon_fbp_hu_f32.raw"), hir=load("recon_hir_hu_f32.raw"),
            shape=shape)
end

@info "Loading recon volumes…"
runs = [load_recon(r) for r in RATES]
shape = runs[1].shape
mid = shape[3] ÷ 2 + 1
@info "  recon shape = $shape, mid-slice z=$mid"

# Two windows: soft-tissue / CTA
windows = [
    (name="soft_tissue", range=(-100.0, 400.0)),  # W=500, L=150
    (name="cta",         range=(-50.0, 700.0)),   # W=750, L=325
]

algos = [("FBP", :fbp), ("HIR", :hir)]
fig = CM.Figure(size=(2000, 700 * length(windows)))
for (row, w) in enumerate(windows)
    col = 1
    for (algo_label, algo_field) in algos
        for run in runs
            img = getfield(run, algo_field)[:, :, mid]
            ax = CM.Axis(fig[row, col]; aspect=CM.DataAspect(),
                title="rate=$(run.rate) — $(algo_label) — $(w.name) " *
                      "(W=$(round(Int, w.range[2] - w.range[1])), L=$(round(Int, (w.range[1] + w.range[2])/2)))",
                titlesize=16)
            CM.heatmap!(ax, img; colormap=:grays, colorrange=w.range)
            CM.hidedecorations!(ax)
            col += 1
        end
    end
    # Colorbar at the end of the row
    CM.Colorbar(fig[row, col]; colormap=:grays, colorrange=w.range,
                label="HU", width=14, labelsize=14)
end

out_path = joinpath(ROOT, "phantom_ct_input", "uci_compare.png")
CM.save(out_path, fig; px_per_unit=2)
@info "Wrote $out_path"

# Also produce a heart-zoom version: crop center 200×200 of mid slice.
center = size(runs[1].fbp, 1) ÷ 2
crop = (center - 80):(center + 80)
fig2 = CM.Figure(size=(2200, 1100))
for (row, w) in enumerate(windows)
    col = 1
    for (algo_label, algo_field) in algos
        for run in runs
            img = getfield(run, algo_field)[crop, crop, mid]
            ax = CM.Axis(fig2[row, col]; aspect=CM.DataAspect(),
                title="rate=$(run.rate) — $(algo_label) — $(w.name) heart zoom",
                titlesize=18)
            CM.heatmap!(ax, img; colormap=:grays, colorrange=w.range)
            CM.hidedecorations!(ax)
            col += 1
        end
    end
    CM.Colorbar(fig2[row, col]; colormap=:grays, colorrange=w.range,
                label="HU", width=18, labelsize=18)
end
out_zoom = joinpath(ROOT, "phantom_ct_input", "uci_compare_heart_zoom.png")
CM.save(out_zoom, fig2; px_per_unit=2)
@info "Wrote $out_zoom"
