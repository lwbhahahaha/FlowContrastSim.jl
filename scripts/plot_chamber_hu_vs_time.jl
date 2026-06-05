#!/usr/bin/env julia
# plot_chamber_hu_vs_time.jl — render HU(t) curves from the time-sweep CSV.
#
# CLI:  julia --project=. plot_chamber_hu_vs_time.jl OUTPUT_DIR
#       (OUTPUT_DIR has chamber_hu_timeseries.csv and bae_v2_traces.csv)

import Pkg
# Lives in FlowContrastSim.jl/scripts/ but uses CairoMakie from the BS scripts
# project (which is also where time_series_dect.jl runs).
Pkg.activate(joinpath(@__DIR__, "..", "..", "phantom_ct_input", "run_basis_sim"))
using CairoMakie
using DelimitedFiles

length(ARGS) >= 1 || error("Usage: plot_chamber_hu_vs_time.jl OUTPUT_DIR")
out = abspath(ARGS[1])
isdir(out) || error("not a dir: $out")
csv = joinpath(out, "chamber_hu_timeseries.csv")
isfile(csv) || error("missing CSV: $csv")

data, header_row = readdlm(csv, ','; header=true)
hdr = vec(String.(header_row))
t = Float64.(data[:, 1])

# Find HU columns
col_of(name) = findfirst(==(name), hdr)
ch_show = [("DA", :DA, :firebrick),
           ("LV", :LV_blood_pool, :forestgreen),
           ("RV", :RV_blood_pool, :royalblue),
           ("PA", :pulm_artery,   :goldenrod),
           ("PV", :pulm_veins,    :purple),
           ("GreatVess.", :great_vessels, :black),]

fig = Figure(size=(1400, 600), fontsize=14)
ax1 = Axis(fig[1,1],
    title="UCI Triphasic 100 kg (66/33/30 @ 5/1.5/2.5)  —  Bae V2 PBPK  —  120 kVp BasisSim",
    xlabel="Time from injection start (s)", ylabel="HU (recon)")
ax2 = Axis(fig[1,2],
    title="Bae V2 predicted [I] (mgI/mL)",
    xlabel="Time from injection start (s)", ylabel="C_iodine (mgI/mL)")

for (lbl, sym, c) in ch_show
    col = col_of("HU_$(String(sym))")
    if col === nothing
        @warn "missing HU column for $sym"; continue
    end
    hu = Float64.(data[:, col])
    lines!(ax1, t, hu, label=lbl, color=c, linewidth=2.5)
    scatter!(ax1, t, hu, color=c, markersize=5)

    col_C = col_of("bae_C_$(String(sym))_mgI_ml")
    if col_C !== nothing
        C = Float64.(data[:, col_C])
        lines!(ax2, t, C, label=lbl, color=c, linewidth=2.5, linestyle=:dash)
    end
end
hlines!(ax1, [40], color=:gray, linestyle=:dash, linewidth=1)
hlines!(ax1, [180, 200], color=:black, linestyle=:dot, linewidth=1)
axislegend(ax1, position=:rt)
axislegend(ax2, position=:rt)
out_png = joinpath(out, "chamber_hu_vs_time.png")
save(out_png, fig)
println("Saved $(out_png)")
