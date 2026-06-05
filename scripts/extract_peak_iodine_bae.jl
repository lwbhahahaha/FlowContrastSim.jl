#!/usr/bin/env julia
#
# extract_peak_iodine_bae.jl — like extract_peak_iodine.jl, but feeds each
# coronary tree's Taylor-Aris PDE the Bae 1998 PBPK aortic-root TAC instead
# of the legacy gamma-variate root input.
#
# Pipeline:
#   1. Run Bae for the UCI triphasic protocol (49/X/30 mL @ 5/Y/2.5 mL/s).
#      Y is the phase-2 rate sweep axis (CLI ARG).
#   2. Take C_aorta(t) as the AIF for all 3 coronary trees.
#   3. For each tree: hemo + simulate_contrast(aif=…) → segment TACs.
#   4. Global peak time = argmax of total iodine mass across trees.
#   5. Per-segment peak iodine + arrival times → {name}_peak_iodine.f32,
#      {name}_arrival_time.f32 (V4 cross-product voxelizer reads these).
#   6. peak_metadata.toml: includes a [chamber_concentrations] section with
#      Bae C_RV/C_LV/C_aorta/C_PA/C_PV at scan time (= t_trigger + scan_delay).
#      The downstream chamber-patch script (add_chambers_to_phantom.jl)
#      reads this section to place iodine in the XCAT chamber labels.
#
# CLI:
#   julia --project=. scripts/extract_peak_iodine_bae.jl \
#         TREE_DIR  OUTPUT_DIR  PHASE2_RATE_ML_S  [SCAN_DELAY_S=6.1]
#
# Default phase 1 = 49 mL @ 5 mL/s, phase 3 = 30 mL saline @ 2.5 mL/s,
# weight 70 kg / 173 cm / male, GE Revolution 120 kVp.

using FlowContrastSim
using Printf
using TOML

# ── Inputs ───────────────────────────────────────────────────────────
const WEIGHT_KG  = 100.0
const HEIGHT_CM  = 173.0
const PHASE1_VOL = 66.0;  const PHASE1_RATE = 5.0
const PHASE2_VOL = 33.0   # fixed; phase2 RATE is the sweep axis
const PHASE3_VOL = 30.0;  const PHASE3_RATE = 2.5
const CONC_MGI_ML= 370.0   # Isovue 370
const KVP        = 120.0
const TRIGGER_HU = 150.0
const BASELINE_HU= 40.0

# ── Hemo + flow config (same physics, but root input becomes Bae AIF) ──
const ROOT_PRESSURE_MMHG    = 100.0
const TERMINAL_PRESSURE_MMHG= 15.0
const HEMATOCRIT            = 0.45
const CAP_BED_R             = 0.12     # hyperemic
const DT                    = 0.1
const T_END                 = 60.0     # cover scan + a bit
const CONTRAST_MIN_DIAM_UM  = 50.0     # sparse threshold
# Mass tables (per-tree territory g) — required for cap-bed scaling.
const TERRITORY_MASSES_G = Dict("LAD" => 58.9, "LCX" => 60.9, "RCA" => 63.8)

# ── BFS helpers (copied from extract_peak_iodine.jl) ─────────────────
function bfs_segment_order(tree)
    n_segs = length(tree.segment_start)
    order = Vector{Int}(undef, 0); sizehint!(order, n_segs)
    queue = Int[]
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && push!(queue, seg)
    end
    visited = falses(n_segs)
    head = 1
    while head <= length(queue)
        s = queue[head]; head += 1
        (s < 1 || s > n_segs || visited[s]) && continue
        visited[s] = true; push!(order, s)
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
    order = bfs_segment_order(tree)
    @inbounds for s in order
        simulated[s] && continue
        start_v = tree.segment_start[s]
        pseg = tree.incoming_segment[start_v]
        pseg > 0 && (peak[s] = peak[pseg])
    end
    return peak
end

# ── Linear interp on saved Bae time grid ──────────────────────────────
function _interp_grid(times, C, t)
    t <= times[1] && return 0.0
    t >= times[end] && return C[end]
    i = searchsortedlast(times, t)
    frac = (t - times[i]) / (times[i+1] - times[i])
    return C[i]*(1-frac) + C[i+1]*frac
end

# ────────────────────────────────────────────────────────────────────
function main()
    if length(ARGS) < 3
        println("Usage: julia --project=. scripts/extract_peak_iodine_bae.jl  TREE_DIR  OUTPUT_DIR  PHASE2_RATE_ML_S  [SCAN_DELAY_S=6.1]")
        exit(1)
    end
    tree_dir  = ARGS[1]
    out_dir   = ARGS[2]
    phase2_rate = parse(Float64, ARGS[3])
    scan_delay  = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 6.1
    isdir(tree_dir) || error("tree_dir not found: $tree_dir")
    mkpath(out_dir)

    @info "Bae UCI triphasic — phase2 rate $(phase2_rate) mL/s, scan_delay $(scan_delay) s after trigger"

    # ── 1. Run Bae central circulation ──
    patient = Patient(weight_kg=WEIGHT_KG, height_cm=HEIGHT_CM)
    protocol = TriphasicProtocol(
        weight_kg=WEIGHT_KG, contrast_concentration_mgI_ml=CONC_MGI_ML,
        phase1_volume_ml=PHASE1_VOL, phase1_rate_ml_s=PHASE1_RATE,
        phase2_volume_ml=PHASE2_VOL, phase2_rate_ml_s=phase2_rate,
        phase3_volume_ml=PHASE3_VOL, phase3_rate_ml_s=PHASE3_RATE,
        phase2_dilution=1.0, phase3_dilution=0.0,
    )
    bae = simulate_central_circulation(patient, protocol; tspan=(0.0, T_END), dt_save=0.05)
    t_trigger = bolus_trigger_time(bae; threshold_delta_HU=TRIGGER_HU, kvp=KVP)
    t_scan    = t_trigger + scan_delay
    @info "  trigger @ $(round(t_trigger,digits=2))s, scan @ $(round(t_scan,digits=2))s"

    # Sample AIF onto the flow simulation's time grid
    times_flow = collect(0.0:DT:T_END)
    aif_vec = [_interp_grid(bae.times, bae.C_aorta, t) for t in times_flow]
    aif_peak = maximum(aif_vec); aif_peak_t = (argmax(aif_vec) - 1) * DT
    @info "  AIF (Bae C_aorta on flow grid): peak $(round(aif_peak,digits=3)) mgI/mL @ t=$(round(aif_peak_t,digits=2))s"

    # ── 2. Discover trees ──
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
    @info "Trees: $(join([n for (n,_) in tree_paths], ", "))"

    # ── 3. Hemo + Bae-AIF contrast per tree ──
    root_p_pa = ROOT_PRESSURE_MMHG * 133.322
    term_p_pa = TERMINAL_PRESSURE_MMHG * 133.322
    trees = Dict{String, Any}()
    hemos = Dict{String, Any}()
    crs   = Dict{String, Any}()

    for (name, path) in tree_paths
        t0 = time()
        @info "[$name] load_tree ..."
        tree = FlowContrastSim.load_tree(name, path)
        @info "[$name] loaded $(length(tree.segment_start)) segs in $(round(time()-t0,digits=1))s"

        mass_g = get(TERRITORY_MASSES_G, name, 0.0)
        t1 = time()
        hemo = FlowContrastSim.compute_hemodynamics(tree;
            root_pressure=root_p_pa, terminal_pressure=term_p_pa,
            hematocrit=HEMATOCRIT,
            capillary_bed_R_per_100g_mmHgmin_ml=CAP_BED_R,
            territory_mass_g=mass_g)
        root_flow = 0.0
        for c in tree.children[tree.root_vertex]
            seg = tree.incoming_segment[c]
            seg != 0 && (root_flow += hemo.segment_flow[seg])
        end
        @info "[$name] hemo $(round(time()-t1,digits=1))s  root_flow=$(round(root_flow*60e6,digits=1)) mL/min"

        t2 = time()
        cr = FlowContrastSim.simulate_contrast(tree, hemo;
            dt=DT, t_end=T_END, aif=aif_vec,
            min_diameter_um=CONTRAST_MIN_DIAM_UM)
        peakC = maximum(cr.concentration)
        @info "[$name] contrast $(round(time()-t2,digits=1))s  peak C=$(round(peakC,digits=3)) mgI/mL"

        trees[name] = tree;  hemos[name] = hemo;  crs[name] = cr
        GC.gc()
    end

    # ── 4. Find global peak time (across trees) ──
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
    @info "global peak @ t=$(round(peak_time,digits=2))s (timestep $(peak_ti)/$(nt))"

    # Note: scan-time chamber values come from Bae directly, not from peak_time
    # of the tree contrast (those may differ — the tree's contrast wave peak is
    # downstream of the aorta peak by the tree's transit time, but the chambers
    # see the BAE-time peak shape directly).
    C_aorta_scan = _interp_grid(bae.times, bae.C_aorta, t_scan)
    C_RV_scan    = _interp_grid(bae.times, bae.C_RV,    t_scan)
    C_LV_scan    = _interp_grid(bae.times, bae.C_LV,    t_scan)
    C_PA_scan    = _interp_grid(bae.times, bae.C_pulm_artery, t_scan)
    C_PV_scan    = _interp_grid(bae.times, bae.C_pulm_vein,   t_scan)
    slope = hu_per_mgI_ml(KVP)
    @info "chamber C @ t_scan (HU at $(KVP) kVp, baseline $(BASELINE_HU)):"
    for (nm, c) in (("aorta",C_aorta_scan), ("RV",C_RV_scan), ("LV",C_LV_scan),
                    ("PA",C_PA_scan), ("PV",C_PV_scan))
        @printf("  %-5s  C=%.2f mgI/mL  HU=%.0f\n", nm, c, BASELINE_HU + slope*c)
    end

    # ── 5. Per-segment peak iodine + arrival files ──
    max_seg_c_global = 0.0
    per_tree_meta = Dict{String, Any}()
    for (name, tree) in trees
        cr = crs[name]
        n_segs = length(tree.segment_start)
        t0 = time()
        peak_iodine = build_peak_iodine_vector(tree, cr, peak_ti)
        max_seg_c = Float64(maximum(peak_iodine))
        max_seg_c_global = max(max_seg_c_global, max_seg_c)
        nonzero = count(>(0f0), peak_iodine)

        out_bin = joinpath(out_dir, "$(lowercase(name))_peak_iodine.f32")
        open(out_bin, "w") do io; write(io, peak_iodine); end

        arrival_f32 = Vector{Float32}(undef, n_segs)
        @inbounds for s in 1:n_segs
            a = cr.arrival_s[s]
            arrival_f32[s] = isfinite(a) ? Float32(a) : Inf32
        end
        out_arr = joinpath(out_dir, "$(lowercase(name))_arrival_time.f32")
        open(out_arr, "w") do io; write(io, arrival_f32); end

        per_tree_meta[name] = (n_segs=n_segs, n_with_iodine=nonzero,
                                max_iodine_mg_per_mL=max_seg_c)
        @info "[$name] $(n_segs) segs, $(nonzero) with iodine>0, max C=$(round(max_seg_c,digits=3)) ($(round(time()-t0,digits=1))s)"
    end

    # ── 6. peak_metadata.toml — bumped iodine_max to cover chambers ──
    chamber_max = max(C_aorta_scan, C_RV_scan, C_LV_scan, C_PA_scan, C_PV_scan)
    iodine_max = max(max_seg_c_global, chamber_max) * 1.05    # 5% headroom
    @info "iodine_max for cross-product encoding = $(round(iodine_max,digits=3)) mgI/mL  (tree max=$(round(max_seg_c_global,digits=3)), chamber max=$(round(chamber_max,digits=3)))"

    meta_path = joinpath(out_dir, "peak_metadata.toml")
    open(meta_path, "w") do io
        println(io, "# peak_metadata.toml — Bae 1998 PBPK AIF, UCI triphasic protocol")
        println(io, "# Generated by FlowContrastSim/scripts/extract_peak_iodine_bae.jl")
        println(io)
        println(io, "[peak]")
        @printf(io, "time_s = %.6f\n", peak_time)
        @printf(io, "timestep = %d\n", peak_ti)
        @printf(io, "n_timesteps = %d\n", nt)
        @printf(io, "total_iodine_mass_mg = %.6f\n", total_mass_mg[peak_ti])
        @printf(io, "max_iodine_concentration_mg_per_mL = %.6f\n", iodine_max)
        @printf(io, "tree_max_iodine_mg_per_mL = %.6f\n", max_seg_c_global)
        @printf(io, "chamber_max_iodine_mg_per_mL = %.6f\n", chamber_max)
        println(io)
        println(io, "[protocol]")
        println(io, "name = \"UCI triphasic\"")
        @printf(io, "weight_kg = %.1f\n", WEIGHT_KG)
        @printf(io, "contrast_concentration_mgI_ml = %.1f\n", CONC_MGI_ML)
        @printf(io, "phase1_volume_ml = %.1f\n", PHASE1_VOL)
        @printf(io, "phase1_rate_ml_s = %.2f\n", PHASE1_RATE)
        @printf(io, "phase2_volume_ml = %.1f\n", PHASE2_VOL)
        @printf(io, "phase2_rate_ml_s = %.2f\n", phase2_rate)
        @printf(io, "phase3_volume_ml = %.1f\n", PHASE3_VOL)
        @printf(io, "phase3_rate_ml_s = %.2f\n", PHASE3_RATE)
        println(io)
        println(io, "[scanner]")
        @printf(io, "kVp = %.1f\n", KVP)
        @printf(io, "trigger_threshold_delta_HU = %.1f\n", TRIGGER_HU)
        @printf(io, "scan_delay_s = %.2f\n", scan_delay)
        @printf(io, "baseline_HU = %.1f\n", BASELINE_HU)
        @printf(io, "trigger_time_s = %.4f\n", t_trigger)
        @printf(io, "scan_time_s = %.4f\n", t_scan)
        println(io)
        println(io, "[chamber_concentrations]")
        println(io, "# Bae 1998 PBPK at scan_time_s. Units: mg I/mL of blood.")
        println(io, "# Consumed by add_chambers_to_phantom.jl to encode")
        println(io, "# cross-product labels for XCAT chamber blood pools.")
        @printf(io, "aorta_root_mgI_ml = %.6f\n", C_aorta_scan)
        @printf(io, "right_heart_mgI_ml = %.6f\n", C_RV_scan)
        @printf(io, "left_heart_mgI_ml = %.6f\n", C_LV_scan)
        @printf(io, "pulm_artery_mgI_ml = %.6f\n", C_PA_scan)
        @printf(io, "pulm_vein_mgI_ml = %.6f\n", C_PV_scan)
        println(io)
        println(io, "[trees]")
        for name in sort(collect(keys(per_tree_meta)))
            m = per_tree_meta[name]
            println(io, "[trees.$name]")
            @printf(io, "n_segs = %d\n", m.n_segs)
            @printf(io, "n_with_iodine_at_peak = %d\n", m.n_with_iodine)
            @printf(io, "max_iodine_mg_per_mL = %.6f\n", m.max_iodine_mg_per_mL)
            @printf(io, "peak_iodine_file = \"%s_peak_iodine.f32\"\n", lowercase(name))
            @printf(io, "arrival_time_file = \"%s_arrival_time.f32\"\n", lowercase(name))
        end
    end
    @info "wrote $(meta_path)"
end

main()
