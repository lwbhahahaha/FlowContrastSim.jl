#!/usr/bin/env julia
#
# build_dynamic_contrast_viewer.jl — generate the dynamic-contrast HTML viewer
# (Plotly + time slider) for the canonical 3-tree max-dilated coronary set.
#
# Lower-level than run_flow_simulation: loads each tree, does hemo + contrast,
# then calls `build_contrast_viewer` directly with **reduced** budgets so the
# Plotly JSON serialization doesn't OOM on 357.8M-segment trees:
#
#   max_segments_per_branch = 800   (vs 6000 default — top-N by diameter)
#   time_stride             = 10    (vs 3 default — keep ~20 frames)
#
# Trees are loaded and processed one at a time, then `GC.gc()` is called
# aggressively to free the full FlowTree before the next is loaded — only the
# sparse ContrastResult (a few-MB) and the small subset of tree segments
# actually needed for the viewer survive across iterations.
#
# Usage:
#   julia --project=. scripts/build_dynamic_contrast_viewer.jl  [TREE_DIR]  [CONFIG]  [OUTPUT_DIR]

using FlowContrastSim
const FCS = FlowContrastSim
using Printf

# ── args ────────────────────────────────────────────────────────────────────
tree_dir    = length(ARGS) >= 1 ? ARGS[1] : "../VascularTreeSim.jl/output"
config_path = length(ARGS) >= 2 ? ARGS[2] : "configs/coronary_hyperemic.toml"
out_dir     = length(ARGS) >= 3 ? ARGS[3] : "../phantom_ct_input/contrast_viewer"

isdir(tree_dir)    || error("tree_dir not found: $tree_dir")
isfile(config_path) || error("config not found: $config_path")
mkpath(out_dir)

println("tree_dir    = $tree_dir")
println("config_path = $config_path")
println("output_dir  = $out_dir")
flush(stdout)

config = FCS.load_flow_config(config_path)
@printf("contrast_amplitude=%.2f mg/mL  dt=%.3fs  t_end=%.1fs (%d timesteps)\n",
        config.contrast_amplitude, config.dt, config.t_end,
        round(Int, config.t_end / config.dt))
flush(stdout)

# ── tree discovery (same regex as run_flow_simulation) ──────────────────────
all_csvs = filter(f -> endswith(lowercase(f), ".csv"), readdir(tree_dir))
vessel_re = r"^([A-Za-z][A-Za-z0-9]*?)(?:_grown)?(?:_segments)?(?:_\d{8}_\d{6})?\.csv$"
aux = ("domain_points.csv", "chambers_points.csv", "pericardium_points.csv",
       "great_vessels_points.csv", "coronary_arteries_points.csv")
tree_paths = Dict{String, String}()
for f in all_csvs
    f in aux && continue
    m = match(vessel_re, f)
    m === nothing && continue
    bname = uppercase(m.captures[1])
    haskey(tree_paths, bname) && occursin(r"_\d{8}_\d{6}\.csv$", f) && continue
    tree_paths[bname] = joinpath(tree_dir, f)
end
branch_names = sort(collect(keys(tree_paths)))
println("Trees: ", branch_names)
flush(stdout)

# ── per-tree run + collect results in memory ────────────────────────────────
root_pressure_pa    = config.root_pressure_mmhg    * 133.322
terminal_pressure_pa = config.terminal_pressure_mmhg * 133.322

trees   = Dict{String, FCS.FlowTree}()
hemos   = Dict{String, FCS.HemodynamicsResult}()
crs     = Dict{String, FCS.ContrastResult}()

for name in branch_names
    path = tree_paths[name]
    t0 = time()
    @info "[$name] load_tree …"
    tree = FCS.load_tree(name, path)
    @info "[$name] loaded $(length(tree.segment_start)) segs in $(round(time()-t0; digits=1))s"
    flush(stdout)

    mass_g = get(config.territory_masses_g, name, 0.0)
    t1 = time()
    hemo = FCS.compute_hemodynamics(tree;
        root_pressure = root_pressure_pa,
        terminal_pressure = terminal_pressure_pa,
        hematocrit = config.discharge_hematocrit,
        capillary_bed_R_per_100g_mmHgmin_ml = config.capillary_bed_R_per_100g_mmHgmin_ml,
        territory_mass_g = mass_g)
    rf = 0.0
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && (rf += hemo.segment_flow[seg])
    end
    @info "[$name] hemo $(round(time()-t1; digits=1))s  root_flow=$(round(rf * 60e6; digits=1)) mL/min"
    flush(stdout)

    t2 = time()
    cr = FCS.simulate_contrast(tree, hemo;
        dt = config.dt, t_end = config.t_end,
        amplitude = config.contrast_amplitude, t0 = config.contrast_t0,
        tmax = config.contrast_tmax, alpha = config.contrast_alpha,
        t_dispersion_s = config.contrast_t_dispersion_s,
        min_diameter_um = config.contrast_min_diameter_um)
    @info "[$name] contrast $(round(time()-t2; digits=1))s  $(size(cr.concentration, 1)) segs sim'd  peak=$(round(maximum(cr.concentration); digits=2)) mg/mL"
    flush(stdout)

    trees[name]   = tree
    hemos[name]   = hemo
    crs[name]     = cr
    GC.gc()
end

# ── build viewer with REDUCED budgets ───────────────────────────────────────
const MAX_SEG_PER_BRANCH = 800
const TIME_STRIDE        = 10

println()
@info "Building viewer (max_segments_per_branch=$MAX_SEG_PER_BRANCH  time_stride=$TIME_STRIDE)…"
flush(stdout)
GC.gc()
t0 = time()
viewer_path = joinpath(out_dir, "contrast_viewer.html")
FCS.build_contrast_viewer(
    viewer_path, trees, hemos, crs;
    time_stride = TIME_STRIDE,
    title = "Dynamic Contrast Transport (coronary trees, t_peak ≈ 10 s)",
    max_segments_per_branch = MAX_SEG_PER_BRANCH,
)
@info "Viewer written in $(round(time()-t0; digits=1))s"
println("Viewer: $viewer_path")
