using Test
using FlowContrastSim
using StaticArrays
using LinearAlgebra

@testset "FlowContrastSim.jl" begin

    @testset "Pries Viscosity Model" begin
        # At large diameters, relative viscosity should approach bulk value (~3.5 cP / 1.2 cP ≈ 2.9)
        eta_large = pries_viscosity_relative(1000.0)
        @test 2.5 < eta_large < 4.0

        # At D ~ 7 μm, near the Fåhræus-Lindqvist minimum — still > 1.0
        eta_small = pries_viscosity_relative(7.0)
        @test eta_small > 1.0

        # At D ~ 3.5 μm (approaching RBC squeeze), viscosity rises steeply
        eta_squeeze = pries_viscosity_relative(3.5)
        @test eta_squeeze > 2.0

        # Below 2.7 μm, effectively infinite
        eta_tiny = pries_viscosity_relative(2.0)
        @test eta_tiny > 1e5

        # Apparent viscosity should be plasma × relative
        mu = apparent_viscosity(100.0)
        @test mu > PLASMA_VISCOSITY_PA_S
        @test mu < 0.01  # should be in reasonable range (Pa·s)

        # Hematocrit dependence: higher Hct → higher viscosity
        eta_low  = pries_viscosity_relative(50.0; hematocrit=0.30)
        eta_high = pries_viscosity_relative(50.0; hematocrit=0.60)
        @test eta_high > eta_low
    end

    @testset "Gamma Variate Input" begin
        times = collect(0.0:0.1:20.0)
        C = gamma_variate_input(times; amplitude=5.0, t0=0.5, tmax=4.0, alpha=3.0)
        @test length(C) == length(times)
        @test C[1] == 0.0            # before t0
        @test maximum(C) > 3.0       # peak should be near amplitude
        @test maximum(C) <= 5.5      # shouldn't exceed amplitude much
        peak_idx = argmax(C)
        @test 3.0 < times[peak_idx] < 5.0  # peak near tmax
    end

    @testset "Synthetic Tree Hemodynamics" begin
        # Build a minimal synthetic FlowTree (Y-bifurcation)
        #   root(1) → v2 → v3 (terminal)
        #                → v4 (terminal)
        vertices = SVector{3,Float64}[
            SVector(0.0, 0.0, 0.0),  # 1: root
            SVector(1.0, 0.0, 0.0),  # 2: bifurcation
            SVector(2.0, 0.5, 0.0),  # 3: terminal left
            SVector(2.0,-0.5, 0.0),  # 4: terminal right
        ]
        parent_vertex = [0, 1, 2, 2]
        children = [Int[2], Int[3, 4], Int[], Int[]]
        incoming_segment = [0, 1, 2, 3]
        segment_start = [1, 2, 2]
        segment_end   = [2, 3, 4]
        # Diameters: parent=0.04cm (400μm), children=0.03cm (300μm)
        segment_diameter_cm = [0.04, 0.03, 0.03]
        segment_label = ["trunk", "left", "right"]

        tree = FlowTree("TestTree", vertices, parent_vertex, children,
                         incoming_segment, segment_start, segment_end,
                         segment_diameter_cm, segment_label, 1)

        # Compute hemodynamics without target flow
        hemo = compute_hemodynamics(tree)
        @test length(hemo.segment_flow) == 3
        @test all(hemo.segment_flow .> 0)  # all segments should have flow

        # Flow conservation: parent flow ≈ sum of children
        @test hemo.segment_flow[1] ≈ hemo.segment_flow[2] + hemo.segment_flow[3] rtol=1e-6

        # Pressure decreases distally
        @test hemo.pressure_proximal[1] > hemo.pressure_distal[1]
        @test hemo.pressure_proximal[2] > hemo.pressure_distal[2]

        # Symmetric children should have equal flow
        @test hemo.segment_flow[2] ≈ hemo.segment_flow[3] rtol=1e-6

        # Transit time should be positive and finite
        @test all(isfinite.(hemo.transit_time_s))
        @test all(hemo.transit_time_s .> 0)

        # Volume should be positive
        @test all(hemo.segment_volume_m3 .> 0)
    end

    @testset "Contrast Transport" begin
        # Same synthetic tree
        vertices = SVector{3,Float64}[
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(2.0, 0.5, 0.0),
            SVector(2.0,-0.5, 0.0),
        ]
        parent_vertex = [0, 1, 2, 2]
        children = [Int[2], Int[3, 4], Int[], Int[]]
        incoming_segment = [0, 1, 2, 3]
        segment_start = [1, 2, 2]
        segment_end   = [2, 3, 4]
        segment_diameter_cm = [0.04, 0.03, 0.03]
        segment_label = ["trunk", "left", "right"]

        tree = FlowTree("TestTree", vertices, parent_vertex, children,
                         incoming_segment, segment_start, segment_end,
                         segment_diameter_cm, segment_label, 1)

        hemo = compute_hemodynamics(tree)
        cr = simulate_contrast(tree, hemo; dt=0.1, t_end=10.0)

        @test length(cr.times) > 50
        @test size(cr.concentration, 1) == 3   # 3 segments
        @test size(cr.concentration, 2) == length(cr.times)

        # Contrast should arrive (nonzero at later times)
        @test maximum(cr.concentration) > 0.1

        # At t=0, no contrast
        @test all(cr.concentration[:, 1] .== 0.0)
    end

    @testset "FlowConfig TOML" begin
        # Write a temp TOML and load it
        toml_path = tempname() * ".toml"
        open(toml_path, "w") do io
            println(io, """
            root_pressure_mmhg = 100.0
            terminal_pressure_mmhg = 15.0
            discharge_hematocrit = 0.45
            dt = 0.1
            t_end = 20.0
            contrast_amplitude = 5.0
            contrast_t0 = 0.5
            contrast_tmax = 4.0
            contrast_alpha = 3.0
            max_arrival_s = 15.0

            [target_flows_ml_min]
            LAD = 242.0
            LCX = 116.0
            RCA = 214.0
            """)
        end

        config = load_flow_config(toml_path)
        @test config.root_pressure_mmhg == 100.0
        @test config.terminal_pressure_mmhg == 15.0
        @test config.target_flows_ml_min["LAD"] == 242.0
        @test config.target_flows_ml_min["RCA"] == 214.0
        @test config.dt == 0.1

        rm(toml_path)
    end

    @testset "Contrast Viewer Generation" begin
        # Minimal tree for viewer test
        vertices = SVector{3,Float64}[
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(2.0, 0.5, 0.0),
        ]
        parent_vertex = [0, 1, 2]
        children = [Int[2], Int[3], Int[]]
        incoming_segment = [0, 1, 2]
        segment_start = [1, 2]
        segment_end   = [2, 3]
        segment_diameter_cm = [0.04, 0.03]
        segment_label = ["trunk", "branch"]

        tree = FlowTree("T", vertices, parent_vertex, children,
                         incoming_segment, segment_start, segment_end,
                         segment_diameter_cm, segment_label, 1)

        trees = Dict("T" => tree)
        hemo = compute_hemodynamics(tree)
        hemo_results = Dict("T" => hemo)
        cr = simulate_contrast(tree, hemo; dt=0.5, t_end=5.0)
        contrast_results = Dict("T" => cr)

        html_path = tempname() * ".html"
        build_contrast_viewer(html_path, trees, hemo_results, contrast_results)
        @test isfile(html_path)
        @test filesize(html_path) > 500

        rm(html_path)
    end

    @testset "Injection Protocol → AIF" begin
        # ── UniphaseNoChaser: defaults and rectangular pulse shape ──
        p = UniphaseNoChaser(weight_kg=70.0)
        @test p.contrast_concentration_mgI_ml == 370.0
        @test p.contrast_volume_per_kg == 0.5
        @test p.injection_rate_ml_s == 5.0

        dt = 0.1
        t_max = 60.0
        times, flux = injection_profile(p; dt=dt, t_max=t_max)
        @test length(times) == length(flux)
        @test times[1] == 0.0
        @test times[end] ≈ t_max

        total_vol = p.weight_kg * p.contrast_volume_per_kg                  # 35 mL
        duration = total_vol / p.injection_rate_ml_s                        # 7 s
        peak_flux = p.injection_rate_ml_s * p.contrast_concentration_mgI_ml # 1850 mgI/s
        @test maximum(flux) ≈ peak_flux
        @test flux[1] ≈ peak_flux       # rectangular: nonzero from t=0
        @test flux[end] == 0.0          # zero after duration
        # Mass conservation: ∫ flux dt = total injected iodine (mg)
        total_iodine_mg = total_vol * p.contrast_concentration_mgI_ml       # 12950 mg
        @test sum(flux) * dt ≈ total_iodine_mg rtol=1e-3

        # ── Central transit kernel: unit area, peak near requested time ──
        physio = PatientPhysiology()
        kernel = central_transit_kernel(times, physio.central_transit_delay_s,
                                        physio.central_transit_dispersion_s)
        @test sum(kernel) * dt ≈ 1.0 rtol=1e-3                  # mass-preserving
        kernel_peak_t = (argmax(kernel) - 1) * dt
        @test abs(kernel_peak_t - physio.central_transit_delay_s) < 1.0  # within 1 s of T

        # Narrower dispersion → narrower (taller) kernel
        kernel_narrow = central_transit_kernel(times,
            physio.central_transit_delay_s, 1.5)
        @test maximum(kernel_narrow) > maximum(kernel)

        # ── End-to-end: protocol_to_aif mass conservation ──
        # Total ∫ AIF dt should equal total_injected_iodine / cardiac_output.
        aif_times, aif = protocol_to_aif(p, physio; dt=dt, t_max=t_max)
        @test length(aif) == length(aif_times)
        @test all(>=(0.0), aif)
        expected_auc = total_iodine_mg / physio.cardiac_output_ml_s
        # 2% tolerance: small tail leaks past t_max for low dispersion / short t_max.
        @test sum(aif) * dt ≈ expected_auc rtol=2e-2

        # AIF peak time should land in [delay - 2s, delay + duration + 2s]:
        # convolution shifts the rectangular pulse so its centroid lines up
        # with the kernel peak; for short injections the peak is near `delay`.
        aif_peak_t = (argmax(aif) - 1) * dt
        @test physio.central_transit_delay_s - 2.0 <= aif_peak_t <= physio.central_transit_delay_s + duration + 2.0
    end

    @testset "FlowConfig with [injection_protocol]" begin
        toml_path = tempname() * ".toml"
        open(toml_path, "w") do io
            println(io, """
            root_pressure_mmhg = 100.0
            terminal_pressure_mmhg = 15.0
            discharge_hematocrit = 0.45
            dt = 0.1
            t_end = 60.0

            [territory_masses_g]
            LAD = 58.9

            [target_flows_ml_min]
            LAD = 200.0

            [injection_protocol]
            type = "UniphaseNoChaser"
            weight_kg = 80.0
            contrast_concentration_mgI_ml = 350.0
            contrast_volume_per_kg = 0.6
            injection_rate_ml_s = 6.0

            [patient_physiology]
            cardiac_output_ml_s = 90.0
            central_transit_delay_s = 14.0
            central_transit_dispersion_s = 4.0
            """)
        end
        cfg = load_flow_config(toml_path)
        @test cfg.injection_protocol isa UniphaseNoChaser
        @test cfg.injection_protocol.weight_kg == 80.0
        @test cfg.injection_protocol.contrast_concentration_mgI_ml == 350.0
        @test cfg.injection_protocol.contrast_volume_per_kg == 0.6
        @test cfg.injection_protocol.injection_rate_ml_s == 6.0
        @test cfg.patient_physiology.cardiac_output_ml_s == 90.0
        @test cfg.patient_physiology.central_transit_delay_s == 14.0
        @test cfg.patient_physiology.central_transit_dispersion_s == 4.0

        # with_protocol override returns a new config with the new protocol.
        cfg2 = with_protocol(cfg; protocol=UniphaseNoChaser(weight_kg=90.0))
        @test cfg2.injection_protocol.weight_kg == 90.0
        # Original config unchanged (immutability)
        @test cfg.injection_protocol.weight_kg == 80.0

        # Drop protocol → fallback to gamma-variate.
        cfg3 = with_protocol(cfg; protocol=nothing)
        @test cfg3.injection_protocol === nothing

        # Missing [injection_protocol] section ⇒ field is nothing, defaults used.
        toml_path2 = tempname() * ".toml"
        open(toml_path2, "w") do io
            println(io, "dt = 0.1\nt_end = 20.0\n")
        end
        cfg4 = load_flow_config(toml_path2)
        @test cfg4.injection_protocol === nothing
        @test cfg4.patient_physiology.cardiac_output_ml_s == 83.0    # default

        rm(toml_path); rm(toml_path2)

        # Unknown protocol type raises.
        toml_path3 = tempname() * ".toml"
        open(toml_path3, "w") do io
            println(io, """
            [injection_protocol]
            type = "MysteryProtocol"
            """)
        end
        @test_throws ErrorException load_flow_config(toml_path3)
        rm(toml_path3)
    end

    @testset "Biphase / chaser protocols" begin
        dt = 0.1
        t_max = 60.0

        # ── UniphaseWithChaser: two phases back to back ──
        u_chaser = UniphaseWithChaser(
            weight_kg=70.0,
            chaser_dilution=0.30,
        )
        times, flux = injection_profile(u_chaser; dt=dt, t_max=t_max)
        # Phase 1: 70 × 0.5 / 5.0 = 7.0 s  at 5 × 370 = 1850 mgI/s
        # Phase 2: 70 × 0.5 / 5.0 = 7.0 s  at 5 × 370 × 0.30 = 555 mgI/s
        @test maximum(flux) ≈ 1850.0
        @test flux[1] ≈ 1850.0                       # phase 1 starts at t=0
        # At t = 8 s (during phase 2): expect 555 mgI/s
        idx8 = round(Int, 8.0/dt) + 1
        @test flux[idx8] ≈ 555.0 atol=1e-6
        # After 14 s: zero
        idx15 = round(Int, 15.0/dt) + 1
        @test flux[idx15] == 0.0
        # Mass conservation: total iodine = (vol_contrast × conc) + (vol_chaser × conc × dilution)
        total_iodine = 70.0 * 0.5 * 370.0 + 70.0 * 0.5 * 370.0 * 0.30
        @test sum(flux) * dt ≈ total_iodine rtol=1e-3

        # ── BiphaseNoChaser: two contrast phases at different rates ──
        b = BiphaseNoChaser(weight_kg=70.0,
                            phase1_volume_per_kg=0.4,
                            phase1_rate_ml_s=6.0,    # 28 mL / 6 = 4.67 s
                            phase2_volume_per_kg=0.4,
                            phase2_rate_ml_s=3.0)    # 28 mL / 3 = 9.33 s
        _, fb = injection_profile(b; dt=dt, t_max=t_max)
        # Phase 1: 6 × 370 = 2220 mgI/s, for ~4.67 s
        @test maximum(fb) ≈ 2220.0
        @test fb[1] ≈ 2220.0
        # Phase 2 starts ~4.67 s, ends ~14.0 s, at 3 × 370 = 1110 mgI/s
        idx_p2 = round(Int, 7.0/dt) + 1            # mid phase 2
        @test fb[idx_p2] ≈ 1110.0 atol=1e-6
        # Mass conservation
        total_b = 70.0 * (0.4 + 0.4) * 370.0
        @test sum(fb) * dt ≈ total_b rtol=1e-3

        # ── BiphaseWithChaser: three back-to-back phases ──
        bc = BiphaseWithChaser(weight_kg=70.0)
        _, fbc = injection_profile(bc; dt=dt, t_max=t_max)
        # Mass: contrast (phase1+phase2) + chaser × dilution
        total_bc = 70.0 * (0.4 + 0.4) * 370.0 + 70.0 * 0.5 * 370.0 * 0.30
        @test sum(fbc) * dt ≈ total_bc rtol=1e-3
        # Three distinct flux levels in order: 2220 (phase1) → 1110 (phase2) → 333 (chaser × 0.3)
        @test maximum(fbc) ≈ 2220.0
        # Find when each phase is active:
        # Phase1 ends at 28/6 ≈ 4.67 s, phase2 ends at 4.67 + 28/3 ≈ 14.0 s,
        # chaser ends at 14.0 + 35/3 ≈ 25.67 s
        @test fbc[round(Int, 2.0/dt)+1] ≈ 2220.0
        @test fbc[round(Int, 8.0/dt)+1] ≈ 1110.0
        @test fbc[round(Int, 20.0/dt)+1] ≈ 0.30 * 3.0 * 370.0  atol=1e-6

        # ── AIF synthesis works for each ──
        physio = PatientPhysiology()
        for p in (u_chaser, b, bc)
            _, aif = protocol_to_aif(p, physio; dt=dt, t_max=t_max)
            @test all(>=(0.0), aif)
            @test maximum(aif) > 0.0
        end
    end

    @testset "FlowConfig parses each protocol" begin
        for (type_str, fields, check) in [
            ("UniphaseNoChaser", "", p -> p isa UniphaseNoChaser),
            ("UniphaseWithChaser",
             "chaser_volume_per_kg = 0.4\nchaser_dilution = 0.20\n",
             p -> p isa UniphaseWithChaser && p.chaser_dilution == 0.20),
            ("BiphaseNoChaser",
             "phase1_volume_per_kg = 0.3\nphase1_rate_ml_s = 7.0\nphase2_volume_per_kg = 0.5\nphase2_rate_ml_s = 2.5\n",
             p -> p isa BiphaseNoChaser && p.phase1_rate_ml_s == 7.0),
            ("BiphaseWithChaser",
             "phase1_volume_per_kg = 0.3\nphase1_rate_ml_s = 7.0\nphase2_volume_per_kg = 0.5\nphase2_rate_ml_s = 2.5\nchaser_dilution = 0.25\n",
             p -> p isa BiphaseWithChaser && p.chaser_dilution == 0.25 && p.phase2_rate_ml_s == 2.5),
        ]
            toml_path = tempname() * ".toml"
            open(toml_path, "w") do io
                println(io, "[injection_protocol]")
                println(io, "type = \"$(type_str)\"")
                println(io, "weight_kg = 75.0")
                println(io, fields)
            end
            cfg = load_flow_config(toml_path)
            @test check(cfg.injection_protocol)
            rm(toml_path)
        end
    end

    @testset "Peak-only mode" begin
        # Build same chain as PDE testset
        vertices = SVector{3,Float64}[
            SVector(0.0,0.0,0.0), SVector(0.5,0.0,0.0),
            SVector(1.0,0.0,0.0), SVector(1.5,0.0,0.0)]
        tree = FlowTree("Chain", vertices, [0,1,2,3], [Int[2],Int[3],Int[4],Int[]],
                        [0,1,2,3], [1,2,3], [2,3,4],
                        [0.4,0.3,0.2], ["s1","s2","s3"], 1)
        hemo = compute_hemodynamics(tree)
        dt = 0.1; t_end = 60.0
        proto = UniphaseNoChaser(weight_kg=70.0)
        physio = PatientPhysiology()
        _, aif = protocol_to_aif(proto, physio; dt=dt, t_max=t_end)

        # Full result first (reference)
        cr_full = simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif)
        @test size(cr_full.concentration) == (3, length(0:dt:t_end))

        # Peak snapshot at a chosen time
        t_peak = 16.0
        cr_peak = simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif, peak_time_s=t_peak)
        @test size(cr_peak.concentration) == (3, 1)
        @test cr_peak.times == [t_peak]

        # Peak-only value should match the full-time result interpolated at t_peak.
        for s in 1:3
            ti = round(Int, t_peak / dt) + 1
            expected = cr_full.concentration[s, ti]
            @test cr_peak.concentration[s, 1] ≈ expected atol=1e-2 rtol=1e-2
        end

        # arrival_variance_s2 is populated and finite on PDE path.
        @test all(isfinite, cr_peak.arrival_variance_s2)
        @test all(>(0.0), cr_peak.arrival_variance_s2)
        @test length(cr_peak.arrival_variance_s2) == length(tree.segment_start)

        # Out-of-range peak_time_s errors.
        @test_throws ErrorException simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif, peak_time_s=-1.0)
        @test_throws ErrorException simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif, peak_time_s=t_end+1.0)

        # Legacy gamma path does not support peak_only.
        @test_throws ErrorException simulate_contrast(tree, hemo; dt=dt, t_end=t_end, peak_time_s=10.0)
    end

    @testset "arrival_variance_s2 in ContrastResult" begin
        vertices = SVector{3,Float64}[
            SVector(0.0,0.0,0.0), SVector(0.5,0.0,0.0),
            SVector(1.0,0.0,0.0), SVector(1.5,0.0,0.0)]
        tree = FlowTree("Chain", vertices, [0,1,2,3], [Int[2],Int[3],Int[4],Int[]],
                        [0,1,2,3], [1,2,3], [2,3,4],
                        [0.4,0.3,0.2], ["s1","s2","s3"], 1)
        hemo = compute_hemodynamics(tree)
        dt = 0.1; t_end = 60.0

        # PDE path: arrival_variance_s2 populated, agrees with _compute_arrival_moments.
        _, aif = protocol_to_aif(UniphaseNoChaser(weight_kg=70.0), PatientPhysiology(); dt=dt, t_max=t_end)
        cr_pde = simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif)
        _, σ²_ref = FlowContrastSim._compute_arrival_moments(tree, hemo)
        @test cr_pde.arrival_variance_s2 ≈ σ²_ref atol=0 rtol=0   # exact float equality

        # Legacy path: arrival_variance_s2 is Inf for every segment.
        cr_leg = simulate_contrast(tree, hemo; dt=dt, t_end=t_end)
        @test all(==(Inf), cr_leg.arrival_variance_s2)
        @test length(cr_leg.arrival_variance_s2) == length(tree.segment_start)
    end

    @testset "Taylor-Aris segment variance" begin
        # Per-segment σ² = 2 D_eff L / v³,   D_eff = D_mol + (Rv)²/(48 D_mol).
        # Pick numerics where the Taylor term dominates so the test value is
        # determined by R, v, L (the empirical leg) not by D_mol.
        R = 1.0e-3            # 1 mm radius
        v = 0.2               # 0.2 m/s
        L = 5.0e-3            # 5 mm
        Dmol = 1.5e-9

        D_eff_expected = Dmol + (R * v)^2 / (48.0 * Dmol)
        σ²_expected    = 2.0 * D_eff_expected * L / v^3

        σ²_got = FlowContrastSim._segment_taylor_variance(R, v, L; D_mol=Dmol)
        @test σ²_got ≈ σ²_expected rtol=1e-12

        # Degenerate inputs return 0 (no dispersion contribution).
        @test FlowContrastSim._segment_taylor_variance(0.0, v, L; D_mol=Dmol) == 0.0
        @test FlowContrastSim._segment_taylor_variance(R, 0.0, L; D_mol=Dmol) == 0.0
        @test FlowContrastSim._segment_taylor_variance(R, v, 0.0; D_mol=Dmol) == 0.0

        # Scaling check: doubling v at fixed (R, L) in the Taylor-dominated
        # regime divides σ² by 2  (Taylor ∝ v², σ² ∝ v² × L / v³ = L/v ⇒ inv).
        σ²_2v = FlowContrastSim._segment_taylor_variance(R, 2v, L; D_mol=Dmol)
        @test σ²_2v ≈ σ²_got / 2.0 rtol=1e-3   # 1% slack for residual D_mol term
    end

    @testset "Arrival moments — additivity along path" begin
        # 3-segment serial chain (no branching): seg 1 → seg 2 → seg 3.
        # σ² at the midpoint of seg 3 should equal
        #   (σ²_full(1) + σ²_full(2)) + σ²_full(3)/2
        # i.e., everything upstream contributes its FULL variance, plus half
        # of the segment containing the sample point.
        vertices = SVector{3,Float64}[
            SVector(0.0, 0.0, 0.0),
            SVector(0.005, 0.0, 0.0),    # 5 mm later
            SVector(0.010, 0.0, 0.0),    # another 5 mm
            SVector(0.015, 0.0, 0.0),    # another 5 mm  (vertices in cm)
        ]
        # Note vertices above are in cm-units (the loader convention) — 5 mm = 0.5 cm
        vertices = SVector{3,Float64}[
            SVector(0.0, 0.0, 0.0),
            SVector(0.5, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(1.5, 0.0, 0.0),
        ]
        parent_vertex    = [0, 1, 2, 3]
        children         = [Int[2], Int[3], Int[4], Int[]]
        incoming_segment = [0, 1, 2, 3]
        segment_start    = [1, 2, 3]
        segment_end      = [2, 3, 4]
        # Trumpet: same diameter, so v adapts to maintain volume flux.
        segment_diameter_cm = [0.2, 0.2, 0.2]   # 2 mm everywhere
        segment_label = ["s1", "s2", "s3"]
        tree = FlowTree("Chain", vertices, parent_vertex, children,
                        incoming_segment, segment_start, segment_end,
                        segment_diameter_cm, segment_label, 1)
        hemo = compute_hemodynamics(tree)

        Dmol = 1.5e-9
        μ, σ² = FlowContrastSim._compute_arrival_moments(tree, hemo; D_mol=Dmol)

        # Compute reference per-segment FULL variance.
        function seg_full_var(s)
            a = tree.vertices[tree.segment_start[s]]
            b = tree.vertices[tree.segment_end[s]]
            L_m = norm(b - a) * 0.01
            R_m = tree.segment_diameter_cm[s] * 0.5 * 0.01
            Q   = abs(hemo.segment_flow[s])
            v   = Q / (π * R_m^2)
            FlowContrastSim._segment_taylor_variance(R_m, v, L_m; D_mol=Dmol)
        end
        sv1 = seg_full_var(1)
        sv2 = seg_full_var(2)
        sv3 = seg_full_var(3)

        # Midpoint of segment 3 sees ALL of seg1, ALL of seg2, HALF of seg3.
        @test σ²[3] ≈ sv1 + sv2 + sv3 / 2.0 rtol=1e-10
        # Midpoint of segment 2 sees ALL of seg1, HALF of seg2.
        @test σ²[2] ≈ sv1 + sv2 / 2.0 rtol=1e-10
        # Midpoint of segment 1 sees HALF of seg1.
        @test σ²[1] ≈ sv1 / 2.0 rtol=1e-10

        # Mean transit time is additive in the same way.
        τ1 = hemo.transit_time_s[1]; τ2 = hemo.transit_time_s[2]; τ3 = hemo.transit_time_s[3]
        @test μ[3] ≈ τ1 + τ2 + τ3 / 2.0 rtol=1e-10
        @test μ[2] ≈ τ1 + τ2 / 2.0 rtol=1e-10
        @test μ[1] ≈ τ1 / 2.0 rtol=1e-10
    end

    @testset "Gaussian convolution sanity" begin
        # ── Mass conservation: a Gaussian with unit-area input keeps unit area. ──
        dt = 0.05
        t_end = 30.0
        times = collect(0.0:dt:t_end)
        # Narrow box around t=10s with area 1
        aif = zeros(length(times))
        i_lo = round(Int, 9.5 / dt) + 1
        i_hi = round(Int, 10.5 / dt)        # 1 s wide
        aif[i_lo:i_hi] .= 1.0               # ∫ AIF dt = 1
        area_in = sum(aif) * dt
        @test area_in ≈ 1.0 rtol=1e-12

        # σ comfortably smaller than the box width and the distance from the
        # box to the time-window edges, so window truncation is negligible.
        F = FlowContrastSim._gaussian_convolve(aif, dt, 0.7)
        @test sum(F) * dt ≈ area_in rtol=1e-3

        # ── Width grows monotonically with σ ──
        F_wide = FlowContrastSim._gaussian_convolve(aif, dt, 2.0)
        # Peak amplitude is lower (mass spread out)
        @test maximum(F_wide) < maximum(F)
        # FWHM (rough): time interval where F > peak/2 is larger for wider σ
        function fwhm(v)
            pk = maximum(v)
            mask = v .>= pk / 2.0
            return sum(mask) * dt
        end
        @test fwhm(F_wide) > fwhm(F)

        # ── σ ≤ dt/2 returns AIF unchanged (delta limit, avoids quantization noise) ──
        F_delta = FlowContrastSim._gaussian_convolve(aif, dt, 0.01)
        @test F_delta == aif
    end

    @testset "PDE simulate_contrast vs legacy" begin
        # Build a small synthetic chain with realistic flow so Taylor is non-trivial.
        vertices = SVector{3,Float64}[
            SVector(0.0, 0.0, 0.0),
            SVector(0.5, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(1.5, 0.0, 0.0),
        ]
        parent_vertex    = [0, 1, 2, 3]
        children         = [Int[2], Int[3], Int[4], Int[]]
        incoming_segment = [0, 1, 2, 3]
        segment_start    = [1, 2, 3]
        segment_end      = [2, 3, 4]
        segment_diameter_cm = [0.4, 0.3, 0.2]
        segment_label = ["s1", "s2", "s3"]
        tree = FlowTree("Chain", vertices, parent_vertex, children,
                        incoming_segment, segment_start, segment_end,
                        segment_diameter_cm, segment_label, 1)
        hemo = compute_hemodynamics(tree)

        dt = 0.1
        t_end = 60.0
        protocol = UniphaseNoChaser(weight_kg=70.0)
        physio = PatientPhysiology()
        _, aif = protocol_to_aif(protocol, physio; dt=dt, t_max=t_end)

        # Run the PDE path
        cr_pde = simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif)
        @test all(cr_pde.concentration[:, 1] .== 0.0)        # zero at t=0
        @test maximum(cr_pde.concentration) > 0.0
        @test maximum(cr_pde.concentration) < maximum(aif)   # passing through ANY dispersion lowers peak

        # ── Mass conservation: ∫C(seg, t)dt ≈ ∫AIF(t)dt for every segment
        # (no leakage of contrast in the model; Gaussian convolution preserves area).
        aif_auc = sum(aif) * dt
        for s in 1:size(cr_pde.concentration, 1)
            seg_auc = sum(cr_pde.concentration[s, :]) * dt
            # 2% tolerance: tail truncation at t_end (deeper segs have peak later).
            @test seg_auc ≈ aif_auc rtol=0.02
        end

        # ── Wide-σ smearing visible: deeper segment has lower peak than shallow.
        @test maximum(cr_pde.concentration[3, :]) < maximum(cr_pde.concentration[1, :])

        # ── σ → 0 limit (zero D_mol → tiny but nonzero Taylor) should give a
        # near-delta convolution, i.e., the result approaches direct AIF
        # interpolation at the arrival time.
        cr_narrow = simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif, D_mol_m2_s=1.0e-15)
        # With D_mol → 0, Taylor term R²v²/(48 D_mol) blows up — variance gets
        # HUGE. So this isn't the right limit. Instead test the *other* direction:
        # huge D_mol shrinks the Taylor term, making σ approach 0.
        cr_big = simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif, D_mol_m2_s=1.0)
        # In this near-delta limit, the segment-1 trace should match
        # _interp_uniform(aif, times, t - arrival[s]) to high accuracy.
        arrival_mean, _ = FlowContrastSim._compute_arrival_moments(tree, hemo; D_mol=1.0)
        times_vec = collect(0.0:dt:t_end)
        a1 = arrival_mean[1]
        for ti in 1:length(times_vec)
            t_local = times_vec[ti] - a1
            expected = t_local < 0.0 ? 0.0 : FlowContrastSim._interp_uniform(aif, times_vec, t_local)
            @test isapprox(cr_big.concentration[1, ti], expected; atol=1e-3, rtol=1e-3)
        end

        # ── PDE result differs from the legacy gamma path (sanity that we
        # actually changed the physics, not just renamed). Use the same AIF
        # for the legacy path is impossible — it doesn't take an AIF — so
        # we run two configs with no AIF (legacy) vs with AIF (PDE) and
        # confirm peaks land at different times because the dispersion model
        # is qualitatively different.
        cr_legacy = simulate_contrast(tree, hemo; dt=dt, t_end=t_end)   # no aif → legacy
        @test cr_legacy.concentration != cr_pde.concentration
    end

    @testset "simulate_contrast consumes AIF" begin
        # Synthetic tree as in earlier testsets
        vertices = SVector{3,Float64}[
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(2.0, 0.5, 0.0),
            SVector(2.0,-0.5, 0.0),
        ]
        parent_vertex = [0, 1, 2, 2]
        children = [Int[2], Int[3, 4], Int[], Int[]]
        incoming_segment = [0, 1, 2, 3]
        segment_start = [1, 2, 2]
        segment_end   = [2, 3, 4]
        segment_diameter_cm = [0.04, 0.03, 0.03]
        segment_label = ["trunk", "left", "right"]
        tree = FlowTree("T", vertices, parent_vertex, children,
                         incoming_segment, segment_start, segment_end,
                         segment_diameter_cm, segment_label, 1)
        hemo = compute_hemodynamics(tree)

        # Build AIF from a uniphase protocol.
        protocol = UniphaseNoChaser(weight_kg=70.0)
        physio = PatientPhysiology()
        dt = 0.1
        t_end = 60.0
        _, aif = protocol_to_aif(protocol, physio; dt=dt, t_max=t_end)

        cr_aif = simulate_contrast(tree, hemo; dt=dt, t_end=t_end, aif=aif)
        @test maximum(cr_aif.concentration) > 0.0
        # AIF-driven peak should land roughly at the AIF peak (transit times
        # in this tiny synthetic tree are sub-second).
        aif_peak_t = (argmax(aif) - 1) * dt
        @test all(cr_aif.concentration[:, 1] .== 0.0)         # no contrast at t=0
        # Each segment's peak time should be close to the AIF peak.
        for s in 1:size(cr_aif.concentration, 1)
            seg_peak_t = (argmax(cr_aif.concentration[s, :]) - 1) * dt
            @test abs(seg_peak_t - aif_peak_t) < 3.0   # ≤3 s slack for tree transit + dispersion
        end

        # Legacy path (no aif) still works.
        cr_legacy = simulate_contrast(tree, hemo; dt=dt, t_end=20.0)
        @test maximum(cr_legacy.concentration) > 0.0
    end

    @testset "VascularTreeSim CSV Contract" begin
        # VascularTreeSim writes parent_segment_id=0 (not -1) for root segments.
        # FlowContrastSim must recognize non-positive values as "no parent".
        csv_path = tempname() * ".csv"
        open(csv_path, "w") do io
            println(io, "branch,segment_id,parent_segment_id,x1_cm,y1_cm,z1_cm,x2_cm,y2_cm,z2_cm,xmid_cm,ymid_cm,zmid_cm,length_mm,diameter_um,label")
            # Root segment with parent_segment_id=0 (VascularTreeSim convention)
            println(io, "LAD,1,0,0.0,0.0,0.0,1.0,0.0,0.0,0.5,0.0,0.0,10.0,400.0,xcat")
            println(io, "LAD,2,1,1.0,0.0,0.0,2.0,0.5,0.0,1.5,0.25,0.0,11.2,300.0,grown")
            println(io, "LAD,3,1,1.0,0.0,0.0,2.0,-0.5,0.0,1.5,-0.25,0.0,11.2,300.0,grown")
        end
        tree = load_tree("LAD", csv_path)
        @test length(tree.segment_start) == 3
        @test length(tree.vertices) == 4     # root, bifurcation, 2 terminals
        # Root vertex must be the one with no incoming segment
        @test tree.incoming_segment[tree.root_vertex] == 0
        @test tree.parent_vertex[tree.root_vertex] == 0
        # Diameters in cm (μm / 1e4)
        @test tree.segment_diameter_cm[1] ≈ 0.04
        @test tree.segment_diameter_cm[2] ≈ 0.03
        # Hemodynamics must run end-to-end on this tree
        hemo = compute_hemodynamics(tree; target_flow_ml_min=10.0)
        @test hemo.segment_flow[1] > 0
        @test hemo.segment_flow[1] ≈ hemo.segment_flow[2] + hemo.segment_flow[3] rtol=1e-6
        rm(csv_path)
    end
end
