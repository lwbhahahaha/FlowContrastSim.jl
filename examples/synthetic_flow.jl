"""
Example: run flow simulation on a synthetic Y-bifurcation tree.

Usage:
    julia --project=.. examples/synthetic_flow.jl

This demonstrates the step-by-step API without needing CSV files.
"""

using FlowContrastSim
using StaticArrays

# ── Step 1: Build a synthetic FlowTree (Y-bifurcation) ──
#
#   root (0,0,0) ─── (2,0,0) ─┬── (4, 1, 0)  terminal
#                               └── (4,-1, 0)  terminal

vertices = SVector{3,Float64}[
    SVector(0.0, 0.0, 0.0),   # 1: root
    SVector(2.0, 0.0, 0.0),   # 2: bifurcation point
    SVector(4.0, 1.0, 0.0),   # 3: left terminal
    SVector(4.0,-1.0, 0.0),   # 4: right terminal
]

tree = FlowTree(
    "YTree",
    vertices,
    [0, 1, 2, 2],                    # parent_vertex
    [Int[2], Int[3, 4], Int[], Int[]], # children
    [0, 1, 2, 3],                     # incoming_segment
    [1, 2, 2],                        # segment_start
    [2, 3, 4],                        # segment_end
    [0.04, 0.025, 0.025],            # segment_diameter_cm (400μm, 250μm)
    ["trunk", "left", "right"],       # segment_label
    1                                  # root_vertex
)

println("Tree: $(length(tree.segment_start)) segments, $(length(tree.vertices)) vertices")

# ── Step 2: Compute hemodynamics ──

hemo = compute_hemodynamics(tree;
    root_pressure=13332.0,        # 100 mmHg
    terminal_pressure=1999.8,     # 15 mmHg
)

println("\nHemodynamics:")
for s in 1:length(tree.segment_start)
    flow_ml_min = hemo.segment_flow[s] * 60e6
    println("  Segment $(s) ($(tree.segment_label[s])): " *
            "flow=$(round(flow_ml_min; digits=3)) mL/min, " *
            "P=$(round(hemo.pressure_proximal[s]/133.322; digits=1))->" *
            "$(round(hemo.pressure_distal[s]/133.322; digits=1)) mmHg, " *
            "transit=$(round(hemo.transit_time_s[s]; digits=4))s")
end

# Flow conservation check
total_child = hemo.segment_flow[2] + hemo.segment_flow[3]
println("\nFlow conservation: parent=$(round(hemo.segment_flow[1]*60e6; digits=3)), " *
        "children sum=$(round(total_child*60e6; digits=3)) mL/min")

# ── Step 3: Simulate contrast transport ──

cr = simulate_contrast(tree, hemo;
    dt=0.1,
    t_end=15.0,
    amplitude=5.0,      # mg/mL iodine
    t0=0.5,             # bolus start
    tmax=4.0,           # bolus peak
    alpha=3.0,          # shape parameter
)

println("\nContrast transport:")
println("  Time points: $(length(cr.times))")
println("  Peak concentration: $(round(maximum(cr.concentration); digits=2)) mg/mL")

for t_check in [2.0, 5.0, 10.0]
    ti = argmin(abs.(cr.times .- t_check))
    conc_at_t = cr.concentration[:, ti]
    println("  t=$(t_check)s: segment concentrations = $(round.(conc_at_t; digits=2))")
end

# ── Step 4: Build viewer ──

outdir = joinpath(@__DIR__, "..", "output", "example_flow")
mkpath(outdir)

trees = Dict("YTree" => tree)
hemo_results = Dict("YTree" => hemo)
contrast_results = Dict("YTree" => cr)

html_path = joinpath(outdir, "index.html")
build_contrast_viewer(html_path, trees, hemo_results, contrast_results;
    title="Y-Bifurcation Flow Example")

println("\nViewer: $(html_path)")
