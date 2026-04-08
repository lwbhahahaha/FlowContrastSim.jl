# ── Top-level orchestration: load trees, compute flow, simulate contrast ──

"""
    run_flow_simulation(tree_csv_dir, config; output_dir="output") -> NamedTuple

Run the full flow simulation pipeline:
1. Auto-discover and load tree CSVs from directory
2. Compute hemodynamics for each tree
3. Simulate contrast transport
4. Build HTML viewer
5. Return results

CSV files are discovered by matching `*_segments.csv` or `*.csv` in the directory.
Branch names are extracted from filenames (e.g., `lad_grown_segments.csv` -> `LAD`).
"""
function run_flow_simulation(tree_csv_dir::String, config::FlowConfig;
                             output_dir::String="output")
    mkpath(output_dir)

    # ── 1. Auto-discover and load tree CSVs ──
    println("[flow] Discovering tree CSVs in: $(tree_csv_dir)")
    csv_files = filter(f -> endswith(f, ".csv"), readdir(tree_csv_dir))
    if isempty(csv_files)
        error("No CSV files found in $(tree_csv_dir)")
    end

    # Extract branch name from filename
    # Handles: lad_grown_segments.csv, LAD.csv, lad_segments.csv, etc.
    function branch_name_from_file(fname::String)
        base = replace(fname, r"\.csv$"i => "")
        # Remove common suffixes
        base = replace(base, r"_grown_segments$"i => "")
        base = replace(base, r"_segments$"i => "")
        base = replace(base, r"_grown$"i => "")
        return uppercase(base)
    end

    tree_paths = Dict{String,String}()
    for f in csv_files
        bname = branch_name_from_file(f)
        tree_paths[bname] = joinpath(tree_csv_dir, f)
        println("[flow]   Found: $(f) -> branch $(bname)")
    end

    trees = load_trees(tree_paths)

    # ── 2. Compute hemodynamics ──
    root_pressure_pa = config.root_pressure_mmhg * 133.322
    terminal_pressure_pa = config.terminal_pressure_mmhg * 133.322

    println("[flow] Computing hemodynamics...")
    hemo_results = Dict{String, HemodynamicsResult}()
    for (name, tree) in trees
        target = get(config.target_flows_ml_min, name, 0.0)
        hemo = compute_hemodynamics(tree;
            root_pressure=root_pressure_pa,
            terminal_pressure=terminal_pressure_pa,
            hematocrit=config.discharge_hematocrit,
            target_flow_ml_min=target)
        hemo_results[name] = hemo

        n_flowing = count(>(1e-18), hemo.segment_flow)
        root_flow = 0.0
        for c in tree.children[tree.root_vertex]
            seg = tree.incoming_segment[c]
            seg != 0 && (root_flow += hemo.segment_flow[seg])
        end
        println("[flow]   $(name): $(n_flowing)/$(length(hemo.segment_flow)) flowing, root_flow=$(round(root_flow*60*1e6, digits=1)) mL/min (target=$(target))")
    end

    # ── 3. Simulate contrast transport ──
    println("[flow] Simulating contrast transport...")
    contrast_results = Dict{String, ContrastResult}()
    for (name, tree) in trees
        cr = simulate_contrast(tree, hemo_results[name];
                               dt=config.dt,
                               t_end=config.t_end,
                               amplitude=config.contrast_amplitude,
                               t0=config.contrast_t0,
                               tmax=config.contrast_tmax,
                               alpha=config.contrast_alpha,
                               max_arrival_s=config.max_arrival_s)
        contrast_results[name] = cr

        cmax = maximum(cr.concentration)
        for t_check in [5.0, 10.0, 15.0]
            t_check > config.t_end && continue
            ti = argmin(abs.(cr.times .- t_check))
            n_conc = count(>(0.01), cr.concentration[:, ti])
            println("[flow]     $(name) @ t=$(t_check)s: $(n_conc) segments with contrast")
        end
    end

    # ── 4. Build HTML viewer ──
    println("[flow] Building contrast viewer...")
    viewer_path = joinpath(output_dir, "contrast_viewer.html")
    build_contrast_viewer(viewer_path, trees, hemo_results, contrast_results;
                          time_stride=3,
                          title="Dynamic Contrast Transport",
                          max_segments_per_branch=6000)
    println("[flow] Viewer written: $(viewer_path)")

    # ── 5. Print summary and return ──
    println("\n" * "="^60)
    println("FLOW SIMULATION SUMMARY")
    println("="^60)
    println("Pressure: $(config.root_pressure_mmhg) mmHg -> $(config.terminal_pressure_mmhg) mmHg")
    println("Hematocrit: $(config.discharge_hematocrit)")
    for name in sort(collect(keys(trees)))
        root_flow = 0.0
        for c in trees[name].children[trees[name].root_vertex]
            seg = trees[name].incoming_segment[c]
            seg != 0 && (root_flow += hemo_results[name].segment_flow[seg])
        end
        target = get(config.target_flows_ml_min, name, 0.0)
        println("  $(name): $(round(root_flow*60*1e6, digits=1)) mL/min (target=$(target))")
    end
    println("Contrast: gamma-variate bolus, $(config.contrast_amplitude) mg/mL")
    println("="^60)

    return (
        trees=trees,
        hemo_results=hemo_results,
        contrast_results=contrast_results,
        viewer_path=viewer_path,
    )
end
