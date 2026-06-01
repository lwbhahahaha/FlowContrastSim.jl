"""
Side-by-side AIF synthesis for all four injection protocols.

Writes `output/aif_protocols.csv` with one column per protocol so they
can be plotted in any tool. The four protocols share `weight_kg = 70.0`,
`contrast_concentration_mgI_ml = 370.0`, and the healthy resting-adult
`PatientPhysiology()` so peak shape differences come purely from the
injection scheme.

Usage:
    julia --project=. examples/compare_protocols_aif.jl
"""

using FlowContrastSim
using Printf

const DT    = 0.1
const T_END = 60.0
const PHYSIO = PatientPhysiology()

# Same weight + iodine concentration; differences are in the phase structure.
protocols = (
    ("UniphaseNoChaser",
     UniphaseNoChaser(weight_kg=70.0)),
    ("UniphaseWithChaser",
     UniphaseWithChaser(weight_kg=70.0)),
    ("BiphaseNoChaser",
     BiphaseNoChaser(weight_kg=70.0)),
    ("BiphaseWithChaser",
     BiphaseWithChaser(weight_kg=70.0)),
)

# Compute all AIFs first.
aifs = Vector{Vector{Float64}}()
times_ref = Float64[]
for (name, p) in protocols
    times, aif = protocol_to_aif(p, PHYSIO; dt=DT, t_max=T_END)
    push!(aifs, aif)
    isempty(times_ref) && (global times_ref = times)
end

# Summary table.
println(rpad("protocol", 22), " | peak (mgI/mL) | peak_t (s) |  AUC (mgI·s/mL) | injected I (mg) | exp AUC")
println("-" ^ 110)
for ((name, p), aif) in zip(protocols, aifs)
    pk = maximum(aif)
    pt = (argmax(aif) - 1) * DT
    auc = sum(aif) * DT
    # Compute expected injected iodine for each scheme.
    inj_I = if p isa UniphaseNoChaser
        p.weight_kg * p.contrast_volume_per_kg * p.contrast_concentration_mgI_ml
    elseif p isa UniphaseWithChaser
        p.weight_kg * p.contrast_volume_per_kg * p.contrast_concentration_mgI_ml +
        p.weight_kg * p.chaser_volume_per_kg   * p.contrast_concentration_mgI_ml * p.chaser_dilution
    elseif p isa BiphaseNoChaser
        p.weight_kg * (p.phase1_volume_per_kg + p.phase2_volume_per_kg) * p.contrast_concentration_mgI_ml
    elseif p isa BiphaseWithChaser
        p.weight_kg * (p.phase1_volume_per_kg + p.phase2_volume_per_kg) * p.contrast_concentration_mgI_ml +
        p.weight_kg * p.chaser_volume_per_kg * p.contrast_concentration_mgI_ml * p.chaser_dilution
    end
    expected_auc = inj_I / PHYSIO.cardiac_output_ml_s
    @printf("%-22s | %12.3f  | %9.2f  | %14.2f | %14.0f | %7.2f\n",
            name, pk, pt, auc, inj_I, expected_auc)
end

# Write CSV with all curves side by side.
out_dir = joinpath(@__DIR__, "..", "output")
mkpath(out_dir)
csv_path = joinpath(out_dir, "aif_protocols.csv")
open(csv_path, "w") do io
    println(io, "t_s," * join((name for (name, _) in protocols), ","))
    for i in eachindex(times_ref)
        vals = (string(round(aif[i], digits=6)) for aif in aifs)
        @printf(io, "%.3f,%s\n", times_ref[i], join(vals, ","))
    end
end
println("\nAIF curves: $(csv_path)")
