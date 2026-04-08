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
end
