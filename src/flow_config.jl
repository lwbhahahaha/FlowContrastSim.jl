# ── Flow simulation configuration from TOML ──

struct FlowConfig
    root_pressure_mmhg::Float64
    terminal_pressure_mmhg::Float64
    discharge_hematocrit::Float64
    target_flows_ml_min::Dict{String,Float64}   # branch_name => target
    contrast_amplitude::Float64
    contrast_t0::Float64
    contrast_tmax::Float64
    contrast_alpha::Float64
    dt::Float64
    t_end::Float64
    max_arrival_s::Float64
    contrast_min_diameter_um::Float64  # 0 = simulate every segment; >0 = skip smaller (saves memory on 8 μm trees)
end

"""
    load_flow_config(path) -> FlowConfig

Load flow simulation configuration from a TOML file.
"""
function load_flow_config(path::String)
    d = TOML.parsefile(path)

    targets = Dict{String,Float64}()
    if haskey(d, "target_flows_ml_min")
        for (k, v) in d["target_flows_ml_min"]
            targets[k] = Float64(v)
        end
    end

    return FlowConfig(
        Float64(get(d, "root_pressure_mmhg", 100.0)),
        Float64(get(d, "terminal_pressure_mmhg", 15.0)),
        Float64(get(d, "discharge_hematocrit", 0.45)),
        targets,
        Float64(get(d, "contrast_amplitude", 5.0)),
        Float64(get(d, "contrast_t0", 0.5)),
        Float64(get(d, "contrast_tmax", 4.0)),
        Float64(get(d, "contrast_alpha", 3.0)),
        Float64(get(d, "dt", 0.1)),
        Float64(get(d, "t_end", 20.0)),
        Float64(get(d, "max_arrival_s", 15.0)),
        Float64(get(d, "contrast_min_diameter_um", 0.0)),
    )
end
