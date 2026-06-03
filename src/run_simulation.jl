# ── Top-level orchestration: load trees, compute flow, simulate contrast ──

"""
    run_flow_simulation(tree_csv_dir, config; output_dir, injection_protocol, patient_physiology)
        -> NamedTuple

Run the full flow simulation pipeline:
1. Auto-discover and load tree CSVs from directory
2. Compute hemodynamics for each tree
3. Synthesize AIF from `config.injection_protocol` + `config.patient_physiology`
   (when set). Per-run overrides take precedence over the TOML defaults:
       run_flow_simulation(dir, cfg; injection_protocol=UniphaseNoChaser(weight_kg=90.0))
4. Simulate contrast transport
5. Build HTML viewer
6. Return results

CSV files are discovered by matching `*_segments.csv` or `*.csv` in the directory.
Branch names are extracted from filenames (e.g., `lad_grown_segments.csv` -> `LAD`).

Pass `injection_protocol=nothing` to force the legacy gamma-variate input
(useful for ablation studies / backward-compatible reproductions).
"""
function run_flow_simulation(tree_csv_dir::String, config::FlowConfig;
                             output_dir::String="output",
                             injection_protocol::Union{Nothing, AbstractInjectionProtocol, Missing}=missing,
                             patient_physiology::Union{PatientPhysiology, Missing}=missing,
                             peak_time_s::Union{Nothing, Float64}=nothing,
                             central_circulation::Union{Nothing, CentralCirculationResult}=nothing)
    mkpath(output_dir)

    # Resolve protocol / physiology: kwarg wins, then config, then default.
    eff_protocol = injection_protocol === missing ? config.injection_protocol : injection_protocol
    eff_physiology = patient_physiology === missing ? config.patient_physiology : patient_physiology

    # Build AIF. Priority:
    #   1. `central_circulation` kwarg — use Bae C_aorta(t) directly (real physics).
    #   2. `injection_protocol` — use protocol_to_aif gamma-kernel path.
    #   3. Otherwise — legacy gamma-variate root input.
    aif_vec = nothing
    if central_circulation !== nothing
        _, aif_vec = aif_from_central(central_circulation; dt=config.dt, t_max=config.t_end)
        println("[flow] Using Bae PBPK AIF from central_circulation kwarg")
        println("[flow]   Patient: $(central_circulation.patient.weight_kg) kg")
        println("[flow]   Protocol: $(typeof(central_circulation.protocol).name.name)")
        flush(stdout)
    elseif eff_protocol !== nothing
        _, aif_vec = protocol_to_aif(eff_protocol, eff_physiology;
                                      dt=config.dt, t_max=config.t_end)
        peak_val = maximum(aif_vec)
        peak_idx = argmax(aif_vec)
        peak_t   = (peak_idx - 1) * config.dt
        auc      = sum(aif_vec) * config.dt
        println("[flow] Protocol: $(typeof(eff_protocol).name.name)")
        println("[flow]   $(eff_protocol)")
        println("[flow]   Physiology: CO=$(eff_physiology.cardiac_output_ml_s) mL/s, transit peak=$(eff_physiology.central_transit_delay_s)s, dispersion=$(eff_physiology.central_transit_dispersion_s)s")
        println("[flow]   AIF: peak=$(round(peak_val, digits=3)) mgI/mL at t=$(round(peak_t, digits=2))s, AUC=$(round(auc, digits=2)) mgI·s/mL")
        flush(stdout)
    else
        println("[flow] No injection_protocol set — falling back to legacy gamma-variate root input.")
        flush(stdout)
    end

    # ── 1. Auto-discover and load tree CSVs ──
    println("[flow] Discovering tree CSVs in: $(tree_csv_dir)")
    all_csvs = filter(f -> endswith(lowercase(f), ".csv"), readdir(tree_csv_dir))
    if isempty(all_csvs)
        error("No CSV files found in $(tree_csv_dir)")
    end

    # Extract branch name and strip optional "_YYYYMMDD_HHMMSS" timestamp.
    # Matches: lad_segments.csv, lad_grown_segments.csv, lad_segments_20260420_151804.csv,
    #          LAD.csv — but skips auxiliary CSVs (domain_points, coronary_arteries_points, ...)
    vessel_re = r"^([A-Za-z][A-Za-z0-9]*?)(?:_grown)?(?:_segments)?(?:_\d{8}_\d{6})?\.csv$"
    aux_files = ("domain_points.csv", "chambers_points.csv", "pericardium_points.csv",
                 "great_vessels_points.csv", "coronary_arteries_points.csv")

    tree_paths = Dict{String,String}()
    for f in all_csvs
        f in aux_files && continue
        m = match(vessel_re, f)
        m === nothing && continue
        bname = uppercase(m.captures[1])
        # When both "lad_segments.csv" (latest) and "lad_segments_<ts>.csv" (snapshot)
        # exist, prefer the non-timestamped latest version.
        if haskey(tree_paths, bname)
            occursin(r"_\d{8}_\d{6}\.csv$", f) && continue
        end
        tree_paths[bname] = joinpath(tree_csv_dir, f)
        println("[flow]   Found: $(f) -> branch $(bname)")
    end
    isempty(tree_paths) && error("No recognizable vessel-tree CSVs in $(tree_csv_dir)")

    # ── 2. Load + compute per tree (sequential to bound memory on 25M-seg 8μm trees) ──
    root_pressure_pa = config.root_pressure_mmhg * 133.322
    terminal_pressure_pa = config.terminal_pressure_mmhg * 133.322

    trees = Dict{String, FlowTree}()
    hemo_results = Dict{String, HemodynamicsResult}()
    contrast_results = Dict{String, ContrastResult}()

    for name in sort(collect(keys(tree_paths)))
        path = tree_paths[name]
        t0 = time()
        println("[flow] [$(name)] loading ..."); flush(stdout)
        tree = load_tree(name, path)
        println("[flow] [$(name)] loaded in $(round(time()-t0,digits=1))s: $(length(tree.segment_start)) segs, $(length(tree.vertices)) verts"); flush(stdout)

        target = get(config.target_flows_ml_min, name, 0.0)
        mass_g = get(config.territory_masses_g, name, 0.0)
        t1 = time()
        # NOTE: target_flow_ml_min is intentionally NOT passed here. The target
        # is kept in the config for display/comparison only — root_flow must
        # emerge from tree topology + Poiseuille + Pries viscosity + physics-
        # based downstream R (`capillary_bed_R_per_100g_mmHgmin_ml`), not from
        # back-solving to match a target.
        hemo = compute_hemodynamics(tree;
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
        n_flowing = count(>(1e-18), hemo.segment_flow)
        println("[flow] [$(name)] hemodynamics in $(round(time()-t1,digits=1))s: $(n_flowing)/$(length(hemo.segment_flow)) flowing, root_flow=$(round(root_flow*60e6,digits=1)) mL/min (target=$(target))"); flush(stdout)

        t2 = time()
        cr = simulate_contrast(tree, hemo;
                               dt=config.dt,
                               t_end=config.t_end,
                               aif=aif_vec,
                               D_mol_m2_s=config.D_mol_m2_s,
                               peak_time_s=peak_time_s,
                               amplitude=config.contrast_amplitude,
                               t0=config.contrast_t0,
                               tmax=config.contrast_tmax,
                               alpha=config.contrast_alpha,
                               t_dispersion_s=config.contrast_t_dispersion_s,
                               min_diameter_um=config.contrast_min_diameter_um)
        nsim = size(cr.concentration, 1)
        println("[flow] [$(name)] contrast in $(round(time()-t2,digits=1))s; simulated $(nsim) segments; peak=$(round(maximum(cr.concentration),digits=2)) mg/mL"); flush(stdout)

        trees[name] = tree
        hemo_results[name] = hemo
        contrast_results[name] = cr
        GC.gc()
    end

    # ── 4. Build HTML viewer ──
    println("[flow] Building contrast viewer..."); flush(stdout)
    viewer_path = joinpath(output_dir, "contrast_viewer.html")
    build_contrast_viewer(viewer_path, trees, hemo_results, contrast_results;
                          time_stride=3,
                          title="Dynamic Contrast Transport",
                          max_segments_per_branch=6000)
    println("[flow] Viewer written: $(viewer_path)"); flush(stdout)

    # ── 5. Print summary and return ──
    println("\n" * "="^60)
    println("FLOW SIMULATION SUMMARY")
    println("="^60)
    println("Pressure: $(config.root_pressure_mmhg) mmHg -> $(config.terminal_pressure_mmhg) mmHg")
    println("Hematocrit: $(config.discharge_hematocrit)")
    println("Capillary bed R: $(config.capillary_bed_R_per_100g_mmHgmin_ml) mmHg·min/mL/100g")
    for name in sort(collect(keys(trees)))
        root_flow = 0.0
        for c in trees[name].children[trees[name].root_vertex]
            seg = trees[name].incoming_segment[c]
            seg != 0 && (root_flow += hemo_results[name].segment_flow[seg])
        end
        target = get(config.target_flows_ml_min, name, 0.0)
        println("  $(name): $(round(root_flow*60*1e6, digits=1)) mL/min (target=$(target))")
    end
    if eff_protocol !== nothing
        println("Contrast input: $(typeof(eff_protocol).name.name) → AIF (gamma-variate central kernel)")
    else
        println("Contrast input: legacy gamma-variate (amplitude=$(config.contrast_amplitude))")
    end
    println("="^60)

    return (
        trees=trees,
        hemo_results=hemo_results,
        contrast_results=contrast_results,
        viewer_path=viewer_path,
        aif=aif_vec,
        protocol=eff_protocol,
        physiology=eff_physiology,
    )
end
