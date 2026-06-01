# ── Flow simulation configuration from TOML ──

struct FlowConfig
    root_pressure_mmhg::Float64
    terminal_pressure_mmhg::Float64
    discharge_hematocrit::Float64
    # Physics-based downstream resistance (replaces target_flow calibration).
    # See `compute_hemodynamics` for interpretation: at max vasodilation this is
    # the true capillary + venous bed R (≈ 0.15); at rest it lumps arteriolar
    # tone + capillary R (≈ 1.0). Same value applied to all trees.
    capillary_bed_R_per_100g_mmHgmin_ml::Float64
    territory_masses_g::Dict{String,Float64}   # branch_name => mass_g
    target_flows_ml_min::Dict{String,Float64}   # branch_name => target (display only)
    contrast_amplitude::Float64
    contrast_t0::Float64
    contrast_tmax::Float64
    contrast_alpha::Float64
    contrast_t_dispersion_s::Float64    # empirical bolus broadening time scale (s); ≈1–3 in vivo coronary
    dt::Float64
    t_end::Float64
    max_arrival_s::Float64               # legacy, accepted but ignored by simulate_contrast
    contrast_min_diameter_um::Float64  # 0 = simulate every segment; >0 = skip smaller (saves memory on 8 μm trees)
    # Iodine molecular diffusivity (m²/s); enters the Taylor-Aris dispersion
    # term in `simulate_contrast` when the AIF path is active. Default
    # 1.5e-9 is mid-range for iohexol/iomeprol in plasma at 37 °C; sensitivity
    # is logarithmic because the Taylor term R²v²/(48 D_mol) and the
    # molecular term D_mol enter with opposite sign of D_mol-dependence.
    D_mol_m2_s::Float64
    # ── Injection protocol + patient physiology (optional). When
    # `injection_protocol` is set, run_flow_simulation builds AIF via
    # protocol_to_aif and feeds it to simulate_contrast — the contrast_*
    # gamma-variate fields above become fallback only. ──
    injection_protocol::Union{Nothing, AbstractInjectionProtocol}
    patient_physiology::PatientPhysiology
end

# ── Protocol parsing from TOML [injection_protocol] table ─────────────────
#
# The `type` key selects which AbstractInjectionProtocol subtype to build.
# All other fields are forwarded as kwargs. Adding a new protocol means
# adding one more branch here + one more struct + injection_profile method
# in protocol.jl — no other code changes.

function _parse_injection_protocol(d::AbstractDict)
    type_str = get(d, "type", "")
    isempty(type_str) && error("[injection_protocol] requires a `type` key, e.g. type = \"UniphaseNoChaser\"")
    haskey(d, "weight_kg") || error("$(type_str) requires weight_kg")
    if type_str == "UniphaseNoChaser"
        return UniphaseNoChaser(
            weight_kg                     = Float64(d["weight_kg"]),
            contrast_concentration_mgI_ml = Float64(get(d, "contrast_concentration_mgI_ml", 370.0)),
            contrast_volume_per_kg        = Float64(get(d, "contrast_volume_per_kg", 0.5)),
            injection_rate_ml_s           = Float64(get(d, "injection_rate_ml_s", 5.0)),
        )
    elseif type_str == "UniphaseWithChaser"
        return UniphaseWithChaser(
            weight_kg                     = Float64(d["weight_kg"]),
            contrast_concentration_mgI_ml = Float64(get(d, "contrast_concentration_mgI_ml", 370.0)),
            contrast_volume_per_kg        = Float64(get(d, "contrast_volume_per_kg", 0.5)),
            injection_rate_ml_s           = Float64(get(d, "injection_rate_ml_s", 5.0)),
            chaser_volume_per_kg          = Float64(get(d, "chaser_volume_per_kg", 0.5)),
            chaser_rate_ml_s              = Float64(get(d, "chaser_rate_ml_s", 5.0)),
            chaser_dilution               = Float64(get(d, "chaser_dilution", 0.30)),
        )
    elseif type_str == "BiphaseNoChaser"
        return BiphaseNoChaser(
            weight_kg                     = Float64(d["weight_kg"]),
            contrast_concentration_mgI_ml = Float64(get(d, "contrast_concentration_mgI_ml", 370.0)),
            phase1_volume_per_kg          = Float64(get(d, "phase1_volume_per_kg", 0.4)),
            phase1_rate_ml_s              = Float64(get(d, "phase1_rate_ml_s", 6.0)),
            phase2_volume_per_kg          = Float64(get(d, "phase2_volume_per_kg", 0.4)),
            phase2_rate_ml_s              = Float64(get(d, "phase2_rate_ml_s", 3.0)),
        )
    elseif type_str == "BiphaseWithChaser"
        return BiphaseWithChaser(
            weight_kg                     = Float64(d["weight_kg"]),
            contrast_concentration_mgI_ml = Float64(get(d, "contrast_concentration_mgI_ml", 370.0)),
            phase1_volume_per_kg          = Float64(get(d, "phase1_volume_per_kg", 0.4)),
            phase1_rate_ml_s              = Float64(get(d, "phase1_rate_ml_s", 6.0)),
            phase2_volume_per_kg          = Float64(get(d, "phase2_volume_per_kg", 0.4)),
            phase2_rate_ml_s              = Float64(get(d, "phase2_rate_ml_s", 3.0)),
            chaser_volume_per_kg          = Float64(get(d, "chaser_volume_per_kg", 0.5)),
            chaser_rate_ml_s              = Float64(get(d, "chaser_rate_ml_s", 3.0)),
            chaser_dilution               = Float64(get(d, "chaser_dilution", 0.30)),
        )
    else
        error("Unknown injection protocol type: \"$(type_str)\". " *
              "Known: UniphaseNoChaser, UniphaseWithChaser, BiphaseNoChaser, BiphaseWithChaser.")
    end
end

function _parse_patient_physiology(d::AbstractDict)
    return PatientPhysiology(
        cardiac_output_ml_s          = Float64(get(d, "cardiac_output_ml_s", 83.0)),
        central_transit_delay_s      = Float64(get(d, "central_transit_delay_s", 12.0)),
        central_transit_dispersion_s = Float64(get(d, "central_transit_dispersion_s", 3.0)),
    )
end

"""
    load_flow_config(path) -> FlowConfig

Load flow simulation configuration from a TOML file.

If the file has an `[injection_protocol]` section, a concrete subtype of
`AbstractInjectionProtocol` is built and stored in `config.injection_protocol`.
Otherwise the field is `nothing` and the legacy gamma-variate parameters
(`contrast_t0`, `contrast_tmax`, `contrast_amplitude`, `contrast_alpha`) are
used as the root input.

The `[patient_physiology]` section is optional; missing fields fall back to
healthy resting-adult defaults (CO 83 mL/s, central transit delay 12 s,
dispersion 3 s).
"""
function load_flow_config(path::String)
    d = TOML.parsefile(path)

    targets = Dict{String,Float64}()
    if haskey(d, "target_flows_ml_min")
        for (k, v) in d["target_flows_ml_min"]
            targets[k] = Float64(v)
        end
    end

    masses = Dict{String,Float64}()
    if haskey(d, "territory_masses_g")
        for (k, v) in d["territory_masses_g"]
            masses[k] = Float64(v)
        end
    end

    protocol = haskey(d, "injection_protocol") ?
        _parse_injection_protocol(d["injection_protocol"]) : nothing
    physiology = haskey(d, "patient_physiology") ?
        _parse_patient_physiology(d["patient_physiology"]) : PatientPhysiology()

    return FlowConfig(
        Float64(get(d, "root_pressure_mmhg", 100.0)),
        Float64(get(d, "terminal_pressure_mmhg", 15.0)),
        Float64(get(d, "discharge_hematocrit", 0.45)),
        Float64(get(d, "capillary_bed_R_per_100g_mmHgmin_ml", 0.0)),
        masses,
        targets,
        Float64(get(d, "contrast_amplitude", 5.0)),
        Float64(get(d, "contrast_t0", 0.5)),
        Float64(get(d, "contrast_tmax", 4.0)),
        Float64(get(d, "contrast_alpha", 3.0)),
        Float64(get(d, "contrast_t_dispersion_s", 3.0)),
        Float64(get(d, "dt", 0.1)),
        Float64(get(d, "t_end", 20.0)),
        Float64(get(d, "max_arrival_s", 15.0)),
        Float64(get(d, "contrast_min_diameter_um", 0.0)),
        Float64(get(d, "D_mol_m2_s", D_MOL_PLASMA_DEFAULT_M2_S)),
        protocol,
        physiology,
    )
end

"""
    with_protocol(config::FlowConfig; protocol=..., physiology=...) -> FlowConfig

Return a new `FlowConfig` with `injection_protocol` and/or
`patient_physiology` replaced. Use this for per-run overrides without
editing the TOML file, e.g.:

    cfg = load_flow_config("configs/coronary_hyperemic_uniphase.toml")
    cfg = with_protocol(cfg; protocol=UniphaseNoChaser(weight_kg=90.0))

`protocol` accepts `nothing` to drop the protocol and fall back to the
legacy gamma-variate input (useful for ablation studies).
"""
function with_protocol(config::FlowConfig;
                       protocol::Union{Nothing, AbstractInjectionProtocol, Missing}=missing,
                       physiology::Union{PatientPhysiology, Missing}=missing)
    new_protocol = protocol === missing ? config.injection_protocol : protocol
    new_physiology = physiology === missing ? config.patient_physiology : physiology
    return FlowConfig(
        config.root_pressure_mmhg,
        config.terminal_pressure_mmhg,
        config.discharge_hematocrit,
        config.capillary_bed_R_per_100g_mmHgmin_ml,
        config.territory_masses_g,
        config.target_flows_ml_min,
        config.contrast_amplitude,
        config.contrast_t0,
        config.contrast_tmax,
        config.contrast_alpha,
        config.contrast_t_dispersion_s,
        config.dt,
        config.t_end,
        config.max_arrival_s,
        config.contrast_min_diameter_um,
        config.D_mol_m2_s,
        new_protocol,
        new_physiology,
    )
end
