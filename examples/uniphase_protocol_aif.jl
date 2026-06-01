"""
Example: synthesize the AIF from a UniphaseNoChaser injection protocol
and print + plot (as CSV) the resulting curve.

Demonstrates two ways to choose the protocol/physiology per run:
  1. Load a TOML config and use its defaults as-is.
  2. Construct/override the protocol or physiology directly in Julia.

Usage:
    julia --project=. examples/uniphase_protocol_aif.jl

The AIF CSV is written to output/uniphase_aif.csv (t_s, AIF_mgI_per_ml).
"""

using FlowContrastSim
using Printf

# ── 1. Load the example TOML config ──
config_path = joinpath(@__DIR__, "..", "configs", "coronary_hyperemic_uniphase.toml")
cfg = load_flow_config(config_path)

println("Loaded config: $(basename(config_path))")
println("  injection protocol: $(typeof(cfg.injection_protocol).name.name)")
println("    weight_kg                     = $(cfg.injection_protocol.weight_kg)")
println("    contrast_concentration_mgI_ml = $(cfg.injection_protocol.contrast_concentration_mgI_ml)")
println("    contrast_volume_per_kg        = $(cfg.injection_protocol.contrast_volume_per_kg)")
println("    injection_rate_ml_s           = $(cfg.injection_protocol.injection_rate_ml_s)")
println("  patient physiology:")
println("    cardiac_output_ml_s           = $(cfg.patient_physiology.cardiac_output_ml_s)")
println("    central_transit_delay_s       = $(cfg.patient_physiology.central_transit_delay_s)")
println("    central_transit_dispersion_s  = $(cfg.patient_physiology.central_transit_dispersion_s)")

# ── 2. Synthesize AIF and inspect its shape ──
times, aif = protocol_to_aif(cfg.injection_protocol, cfg.patient_physiology;
                              dt=cfg.dt, t_max=cfg.t_end)
peak_t  = (argmax(aif) - 1) * cfg.dt
peak_v  = maximum(aif)
auc     = sum(aif) * cfg.dt
total_I = cfg.injection_protocol.weight_kg *
          cfg.injection_protocol.contrast_volume_per_kg *
          cfg.injection_protocol.contrast_concentration_mgI_ml
expected_auc = total_I / cfg.patient_physiology.cardiac_output_ml_s

println("\nAIF stats")
@printf("  peak:               %.3f mgI/mL at t = %.2f s\n", peak_v, peak_t)
@printf("  AUC:                %.2f mgI·s/mL\n", auc)
@printf("  total injected I:   %.0f mg\n", total_I)
@printf("  expected AUC:       %.2f (= total_I / CO)\n", expected_auc)
@printf("  mass conservation:  %.4f×\n", auc / expected_auc)

# ── 3. Per-run override: change protocol / physiology without editing TOML ──
println("\nOverride example — bump weight 70→90 kg and concentration 370→320:")
custom_protocol = UniphaseNoChaser(
    weight_kg                     = 90.0,
    contrast_concentration_mgI_ml = 320.0,
    contrast_volume_per_kg        = 0.6,
    injection_rate_ml_s           = 6.0,
)
custom_physio = PatientPhysiology(
    cardiac_output_ml_s           = 100.0,    # higher CO, e.g. younger patient
    central_transit_delay_s       = 10.0,
    central_transit_dispersion_s  = 2.5,
)
_, aif2 = protocol_to_aif(custom_protocol, custom_physio;
                          dt=cfg.dt, t_max=cfg.t_end)
@printf("  override peak:    %.3f mgI/mL at t = %.2f s\n",
        maximum(aif2), (argmax(aif2) - 1) * cfg.dt)

# Or: keep TOML defaults but swap just the protocol via `with_protocol`.
cfg_modified = with_protocol(cfg; protocol=custom_protocol)
println("  with_protocol(): weight_kg now $(cfg_modified.injection_protocol.weight_kg)" *
        " (original config unchanged: $(cfg.injection_protocol.weight_kg))")

# ── 4. Dump AIF to CSV for downstream tools / plotting ──
out_dir = joinpath(@__DIR__, "..", "output")
mkpath(out_dir)
csv_path = joinpath(out_dir, "uniphase_aif.csv")
open(csv_path, "w") do io
    println(io, "t_s,aif_mgI_per_ml,aif_override_mgI_per_ml")
    for i in eachindex(times)
        @printf(io, "%.3f,%.6f,%.6f\n", times[i], aif[i], aif2[i])
    end
end
println("\nAIF curves written to: $(csv_path)")
