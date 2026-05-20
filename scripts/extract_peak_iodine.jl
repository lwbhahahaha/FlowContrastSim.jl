#!/usr/bin/env julia
#
# extract_peak_iodine.jl — run hemo + contrast for each tree in TREE_DIR,
# find the global peak time (argmax of total iodine mass summed across trees),
# and write a per-segment iodine concentration vector at that time per tree.
#
# The downstream voxelizer (apply_contrast_at_peak.jl in VascularTreeSim.jl)
# consumes these to build a (blood %, iodine %) cross-product UInt16 phantom.
#
# Per-segment iodine handling:
#   - FlowContrastSim only simulates segments with diameter ≥ config.contrast_min_diameter_um
#     (50 μm by default) — sub-threshold segments aren't in the concentration matrix.
#   - We propagate iodine values down the topology (BFS, parent → child) so every
#     reachable segment gets a value. Capillaries inherit from their parent
#     ≥-threshold segment (plug-flow approximation; reasonable since dispersion
#     dominates over capillary transit time).
#   - Sub-threshold segments that have no simulated ancestor get 0.
#
# Output (in OUTPUT_DIR):
#   {lad,lcx,rca}_peak_iodine.f32   — Float32 binary, one value per segment_id
#                                     indexed by tree.segment_start[1..n_segs]
#                                     (i.e., `arr[s] = peak C in seg s`).
#                                     Size: n_segs × 4 bytes.
#   {lad,lcx,rca}_arrival_time.f32  — Float32 binary, per-segment bolus arrival
#                                     time in seconds (Inf32 = unreachable /
#                                     zero-flow). Same indexing convention.
#                                     Consumed by simple_dynamic_viewer.jl.
#   peak_metadata.toml              — peak_time_s, max_concentration, n_segs per tree
#
# Usage:
#   julia --project=. scripts/extract_peak_iodine.jl  TREE_DIR  CONFIG_TOML  OUTPUT_DIR

using FlowContrastSim
using Printf
using TOML

function bfs_segment_order(tree)
    # BFS in segment-space: parent segment before child segments.
    n_segs = length(tree.segment_start)
    order = Vector{Int}(undef, 0)
    sizehint!(order, n_segs)
    queue = Int[]
    # Seed with root segments
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && push!(queue, seg)
    end
    visited = falses(n_segs)
    head = 1
    while head <= length(queue)
        s = queue[head]; head += 1
        (s < 1 || s > n_segs || visited[s]) && continue
        visited[s] = true
        push!(order, s)
        end_v = tree.segment_end[s]
        for c in tree.children[end_v]
            cseg = tree.incoming_segment[c]
            cseg != 0 && !visited[cseg] && push!(queue, cseg)
        end
    end
    order
end

function build_peak_iodine_vector(tree, cr::FlowContrastSim.ContrastResult, peak_ti::Int)
    n_segs = length(tree.segment_start)
    peak = zeros(Float32, n_segs)
    simulated = falses(n_segs)
    seg_ids = isempty(cr.segment_ids) ? (1:size(cr.concentration, 1)) : cr.segment_ids
    @inbounds for (row, s) in enumerate(seg_ids)
        if 1 <= s <= n_segs
            peak[s] = Float32(cr.concentration[row, peak_ti])
            simulated[s] = true
        end
    end

    # Propagate: non-simulated segments inherit from parent (BFS order ensures
    # the parent is processed first).
    order = bfs_segment_order(tree)
    @inbounds for s in order
        simulated[s] && continue
        start_v = tree.segment_start[s]
        pseg = tree.incoming_segment[start_v]
        if pseg > 0
            peak[s] = peak[pseg]   # inherit (may be 0 if parent also non-sim'd)
        end
    end
    return peak
end

function main()
    if length(ARGS) < 3
        println("Usage: julia --project=. scripts/extract_peak_iodine.jl  TREE_DIR  CONFIG_TOML  OUTPUT_DIR")
        exit(1)
    end
    tree_dir, config_path, out_dir = ARGS[1], ARGS[2], ARGS[3]
    isdir(tree_dir) || error("tree_dir not found: $tree_dir")
    isfile(config_path) || error("config not found: $config_path")
    mkpath(out_dir)

    config = FlowContrastSim.load_flow_config(config_path)
    @printf("[cfg] contrast_amplitude=%.2f mg/mL  t0=%.2f  tmax=%.2f  alpha=%.2f  dt=%.3f  t_end=%.1f\n",
            config.contrast_amplitude, config.contrast_t0, config.contrast_tmax,
            config.contrast_alpha, config.dt, config.t_end)
    @printf("[cfg] contrast_min_diameter_um=%.1f  t_dispersion_s=%.2f\n",
            config.contrast_min_diameter_um, config.contrast_t_dispersion_s)
    flush(stdout)

    # ── 1. Discover trees ──
    all_csvs = filter(f -> endswith(lowercase(f), "_segments.csv"), readdir(tree_dir; join=true))
    aux = ("domain_points.csv", "chambers_points.csv", "pericardium_points.csv",
           "great_vessels_points.csv", "coronary_arteries_points.csv")
    tree_paths = Pair{String,String}[]
    for f in sort(all_csvs)
        basename(f) in aux && continue
        m = match(r"^([A-Za-z][A-Za-z0-9]*?)_segments\.csv$", basename(f))
        m === nothing && continue
        push!(tree_paths, uppercase(m.captures[1]) => f)
    end
    isempty(tree_paths) && error("No *_segments.csv found in $tree_dir")
    @printf("[trees] %d:  %s\n", length(tree_paths),
            join([n for (n, _) in tree_paths], ", "))
    flush(stdout)

    # ── 2. Run hemo + contrast for each tree; keep results in memory ──
    root_pressure_pa = config.root_pressure_mmhg * 133.322
    terminal_pressure_pa = config.terminal_pressure_mmhg * 133.322

    trees   = Dict{String, Any}()
    hemos   = Dict{String, Any}()
    crs     = Dict{String, Any}()

    for (name, path) in tree_paths
        t0 = time()
        @printf("\n[%s] load_tree …\n", name); flush(stdout)
        tree = FlowContrastSim.load_tree(name, path)
        @printf("[%s] loaded %d segs in %.1fs\n", name, length(tree.segment_start), time()-t0); flush(stdout)

        mass_g = get(config.territory_masses_g, name, 0.0)
        t1 = time()
        hemo = FlowContrastSim.compute_hemodynamics(tree;
            root_pressure=root_pressure_pa,
            terminal_pressure=terminal_pressure_pa,
            hematocrit=config.discharge_hematocrit,
            capillary_bed_R_per_100g_mmHgmin_ml=config.capillary_bed_R_per_100g_mmHgmin_ml,
            territory_mass_g=mass_g)
        root_flow = 0.0
        for c in tree.children[tree.root_vertex]
            seg = tree.incoming_segment[c]
            seg != 0 && (root_flow += hemo.segment_flow[seg])
        end
        @printf("[%s] hemo %.1fs  root_flow=%.1f mL/min\n",
                name, time()-t1, root_flow*60e6); flush(stdout)

        t2 = time()
        cr = FlowContrastSim.simulate_contrast(tree, hemo;
            dt=config.dt, t_end=config.t_end,
            amplitude=config.contrast_amplitude,
            t0=config.contrast_t0, tmax=config.contrast_tmax, alpha=config.contrast_alpha,
            t_dispersion_s=config.contrast_t_dispersion_s,
            min_diameter_um=config.contrast_min_diameter_um)
        nsim = size(cr.concentration, 1)
        peakC = maximum(cr.concentration)
        finite_arr = filter(isfinite, cr.arrival_s)
        if !isempty(finite_arr)
            arr_min, arr_med, arr_max = minimum(finite_arr),
                                        finite_arr[clamp(div(length(finite_arr), 2), 1, length(finite_arr))],
                                        maximum(finite_arr)
            @printf("[%s] contrast %.1fs  simulated %d segs (≥%.1f μm)  peak C=%.2f mg/mL  arrival: min=%.2fs median=%.2fs max=%.2fs (%d/%d reachable)\n",
                    name, time()-t2, nsim, config.contrast_min_diameter_um, peakC,
                    arr_min, arr_med, arr_max, length(finite_arr), length(cr.arrival_s)); flush(stdout)
        else
            @printf("[%s] contrast %.1fs  simulated %d segs  peak C=%.2f mg/mL  NO REACHABLE SEGMENTS\n",
                    name, time()-t2, nsim, peakC); flush(stdout)
        end

        trees[name] = tree
        hemos[name] = hemo
        crs[name]   = cr
        GC.gc()
    end

    # ── 3. Find global peak time (argmax of total iodine mass across trees) ──
    any_cr = first(values(crs))
    times = any_cr.times
    nt = length(times)
    total_mass_mg = zeros(Float64, nt)
    for (name, cr) in crs
        hemo = hemos[name]
        seg_ids = isempty(cr.segment_ids) ? (1:size(cr.concentration, 1)) : cr.segment_ids
        @inbounds for (row, s) in enumerate(seg_ids)
            vol_mL = hemo.segment_volume_m3[s] * 1e6
            @simd for ti in 1:nt
                total_mass_mg[ti] += cr.concentration[row, ti] * vol_mL
            end
        end
    end
    peak_ti = argmax(total_mass_mg)
    peak_time = times[peak_ti]
    @printf("\n[peak] global peak at t=%.2f s (timestep %d/%d)  total iodine mass=%.3f mg\n",
            peak_time, peak_ti, nt, total_mass_mg[peak_ti])
    flush(stdout)

    # ── 4. Build per-segment peak iodine + write arrival time per tree ──
    max_c_at_peak = 0.0
    per_tree_meta = Dict{String, Any}()
    for (name, tree) in trees
        cr = crs[name]
        n_segs = length(tree.segment_start)
        t0 = time()
        peak_iodine = build_peak_iodine_vector(tree, cr, peak_ti)
        nonzero_after = count(>(0f0), peak_iodine)
        max_seg_c = Float64(maximum(peak_iodine))
        max_c_at_peak = max(max_c_at_peak, max_seg_c)

        out_bin = joinpath(out_dir, "$(lowercase(name))_peak_iodine.f32")
        open(out_bin, "w") do io;  write(io, peak_iodine);  end

        # Arrival time per segment_id (Float32, Inf for unreachable). Length matches
        # full tree so consumers (viewer, voxelizer) can index by segment_id directly.
        arrival_f32 = Vector{Float32}(undef, n_segs)
        @inbounds for s in 1:n_segs
            a = cr.arrival_s[s]
            arrival_f32[s] = isfinite(a) ? Float32(a) : Inf32
        end
        out_arr = joinpath(out_dir, "$(lowercase(name))_arrival_time.f32")
        open(out_arr, "w") do io;  write(io, arrival_f32);  end
        n_reachable = count(isfinite, arrival_f32)

        @printf("[%s] peak_iodine.f32: %d segs (%d with C>0)  max C=%.3f mg/mL  | arrival_time.f32: %d/%d reachable  | %.1fs\n",
                name, n_segs, nonzero_after, max_seg_c, n_reachable, n_segs, time()-t0); flush(stdout)

        per_tree_meta[name] = (
            n_segs        = n_segs,
            n_simulated   = isempty(cr.segment_ids) ? n_segs : length(cr.segment_ids),
            n_with_iodine = nonzero_after,
            n_reachable   = n_reachable,
            max_iodine_mg_per_mL = max_seg_c,
        )
    end

    # ── 5. Write peak_metadata.toml ──
    meta_path = joinpath(out_dir, "peak_metadata.toml")
    open(meta_path, "w") do io
        println(io, "# peak_metadata.toml — output of FlowContrastSim/scripts/extract_peak_iodine.jl")
        println(io)
        println(io, "[peak]")
        @printf(io, "time_s = %.6f\n", peak_time)
        @printf(io, "timestep = %d\n", peak_ti)
        @printf(io, "n_timesteps = %d\n", nt)
        @printf(io, "total_iodine_mass_mg = %.6f\n", total_mass_mg[peak_ti])
        @printf(io, "max_iodine_concentration_mg_per_mL = %.6f\n", max_c_at_peak)
        println(io)
        println(io, "[contrast_bolus]")
        @printf(io, "amplitude_mg_per_mL = %.4f\n", config.contrast_amplitude)
        @printf(io, "t0_s   = %.4f\n", config.contrast_t0)
        @printf(io, "tmax_s = %.4f\n", config.contrast_tmax)
        @printf(io, "alpha  = %.4f\n", config.contrast_alpha)
        @printf(io, "dt_s   = %.4f\n", config.dt)
        @printf(io, "t_end_s= %.4f\n", config.t_end)
        @printf(io, "t_dispersion_s = %.4f\n", config.contrast_t_dispersion_s)
        @printf(io, "contrast_min_diameter_um_simulated = %.4f\n", config.contrast_min_diameter_um)
        println(io)
        println(io, "[trees]")
        for name in sort(collect(keys(per_tree_meta)))
            m = per_tree_meta[name]
            println(io, "[trees.$name]")
            @printf(io, "n_segs                = %d\n", m.n_segs)
            @printf(io, "n_simulated           = %d\n", m.n_simulated)
            @printf(io, "n_reachable           = %d\n", m.n_reachable)
            @printf(io, "n_with_iodine_at_peak = %d\n", m.n_with_iodine)
            @printf(io, "max_iodine_mg_per_mL  = %.6f\n", m.max_iodine_mg_per_mL)
            @printf(io, "peak_iodine_file      = \"%s_peak_iodine.f32\"\n", lowercase(name))
            @printf(io, "arrival_time_file     = \"%s_arrival_time.f32\"\n", lowercase(name))
        end
    end
    @printf("\n[write] %s\n", meta_path)
end

main()
