# ── Contrast injection protocols → AIF synthesis ──
#
# A protocol describes how the contrast bolus is administered into the
# peripheral vein (single-phase, biphase, with/without saline chaser, …).
# Combined with patient physiology (cardiac output + a phenomenological
# central-circulation transit kernel) it produces the AIF (arterial input
# function) at the aorta root — the boundary condition for the in-tree
# Taylor-Aris PDE solver.
#
# Pipeline:
#   protocol            → injection_profile(p)   = mgI/s entering vein
#   ⊛ central_transit_kernel(physiology)         = mgI/s arriving at aorta
#   ÷ cardiac_output_ml_s                        = mgI/mL at aorta = AIF(t)
#
# Adding a new protocol scheme means defining a new <:AbstractInjectionProtocol
# struct, listing its phases as `(volume_ml, rate_ml_s, concentration_mgI_ml)`
# tuples, and adding one elseif to `_parse_injection_protocol` in flow_config.jl.
# `protocol_to_aif` / `central_transit_kernel` are protocol-agnostic.

"""
    AbstractInjectionProtocol

Base type. A subtype encodes one injection scheme and provides an
`injection_profile(p; dt, t_max) -> (times, mgI_per_s)` method.

Currently supported subtypes (all share `weight_kg` + `contrast_concentration_mgI_ml`):

| protocol              | phases                                                    |
|-----------------------|-----------------------------------------------------------|
| `UniphaseNoChaser`    | one rectangular pulse                                      |
| `UniphaseWithChaser`  | contrast pulse → diluted-or-saline chaser                 |
| `BiphaseNoChaser`     | two contrast pulses (typically fast → slow rate)          |
| `BiphaseWithChaser`   | two contrast pulses → chaser                               |
"""
abstract type AbstractInjectionProtocol end

"""
    UniphaseNoChaser(; weight_kg, contrast_concentration_mgI_ml=370.0,
                        contrast_volume_per_kg=0.5, injection_rate_ml_s=5.0)

Single-phase contrast injection, no saline chaser. The injection profile is
a rectangular pulse of duration `weight_kg × contrast_volume_per_kg /
injection_rate_ml_s` at flux `injection_rate_ml_s × contrast_concentration_mgI_ml`.
"""
Base.@kwdef struct UniphaseNoChaser <: AbstractInjectionProtocol
    weight_kg::Float64
    contrast_concentration_mgI_ml::Float64 = 370.0
    contrast_volume_per_kg::Float64 = 0.5
    injection_rate_ml_s::Float64 = 5.0
end

"""
    UniphaseWithChaser(; weight_kg, contrast_concentration_mgI_ml=370.0,
                          contrast_volume_per_kg=0.5, injection_rate_ml_s=5.0,
                          chaser_volume_per_kg=0.5, chaser_rate_ml_s=5.0,
                          chaser_dilution=0.30)

Single-phase contrast injection followed by a saline (or diluted contrast)
chaser. `chaser_dilution` is the iodine fraction in the chaser
(0.0 = pure saline; 0.30 = 30 % contrast / 70 % saline mixed bolus;
1.0 = pure contrast, equivalent to biphase same-rate).

Phases:
1. contrast: `weight_kg × contrast_volume_per_kg` mL at
   `injection_rate_ml_s` mL/s, concentration `contrast_concentration_mgI_ml`.
2. chaser:   `weight_kg × chaser_volume_per_kg` mL at
   `chaser_rate_ml_s` mL/s, concentration
   `contrast_concentration_mgI_ml × chaser_dilution`.
"""
Base.@kwdef struct UniphaseWithChaser <: AbstractInjectionProtocol
    weight_kg::Float64
    contrast_concentration_mgI_ml::Float64 = 370.0
    contrast_volume_per_kg::Float64 = 0.5
    injection_rate_ml_s::Float64 = 5.0
    chaser_volume_per_kg::Float64 = 0.5
    chaser_rate_ml_s::Float64 = 5.0
    chaser_dilution::Float64 = 0.30
end

"""
    BiphaseNoChaser(; weight_kg, contrast_concentration_mgI_ml=370.0,
                       phase1_volume_per_kg=0.4, phase1_rate_ml_s=6.0,
                       phase2_volume_per_kg=0.4, phase2_rate_ml_s=3.0)

Two-phase contrast injection (typical: fast loading then slower sustaining
phase) with no chaser. Both phases use the same contrast concentration.
"""
Base.@kwdef struct BiphaseNoChaser <: AbstractInjectionProtocol
    weight_kg::Float64
    contrast_concentration_mgI_ml::Float64 = 370.0
    phase1_volume_per_kg::Float64 = 0.4
    phase1_rate_ml_s::Float64 = 6.0
    phase2_volume_per_kg::Float64 = 0.4
    phase2_rate_ml_s::Float64 = 3.0
end

"""
    BiphaseWithChaser(; weight_kg, contrast_concentration_mgI_ml=370.0,
                         phase1_volume_per_kg=0.4, phase1_rate_ml_s=6.0,
                         phase2_volume_per_kg=0.4, phase2_rate_ml_s=3.0,
                         chaser_volume_per_kg=0.5, chaser_rate_ml_s=3.0,
                         chaser_dilution=0.30)

Two-phase contrast + chaser. Combines `BiphaseNoChaser` (phases 1+2) with
the chaser semantics of `UniphaseWithChaser`. `chaser_rate_ml_s` is
typically equal to `phase2_rate_ml_s` (the injector continues running) but
is independent here for full flexibility.
"""
Base.@kwdef struct BiphaseWithChaser <: AbstractInjectionProtocol
    weight_kg::Float64
    contrast_concentration_mgI_ml::Float64 = 370.0
    phase1_volume_per_kg::Float64 = 0.4
    phase1_rate_ml_s::Float64 = 6.0
    phase2_volume_per_kg::Float64 = 0.4
    phase2_rate_ml_s::Float64 = 3.0
    chaser_volume_per_kg::Float64 = 0.5
    chaser_rate_ml_s::Float64 = 3.0
    chaser_dilution::Float64 = 0.30
end

"""
    PatientPhysiology(; cardiac_output_ml_s=83.0,
                         central_transit_delay_s=12.0,
                         central_transit_dispersion_s=3.0)

Patient-level parameters that shape the AIF independent of protocol.

* `cardiac_output_ml_s`: cardiac output in mL/s (default 83.0 = 5 L/min).
* `central_transit_delay_s`: peak time of the RV→lung→LV→aorta kernel.
* `central_transit_dispersion_s`: width parameter; smaller = sharper peak.

The kernel is a lumped phenomenological gamma-variate — not a chamber-by-
chamber mechanistic model. We do not solve PDEs in the chambers.
"""
Base.@kwdef struct PatientPhysiology
    cardiac_output_ml_s::Float64 = 83.0
    central_transit_delay_s::Float64 = 12.0
    central_transit_dispersion_s::Float64 = 3.0
end

# ── Common helper: serial rectangular phases → mgI/s vs time ───────────────
#
# Each phase is a NamedTuple `(volume_ml, rate_ml_s, concentration_mgI_ml)`.
# Phases run back-to-back starting at t=0 with no gaps. The flux during
# phase k is `rate_k × concentration_k` (mgI/s), held constant for
# `volume_k / rate_k` seconds.
#
# At phase boundaries we use proportional blending: `flux[i]` is the
# time-average of the instantaneous mass-flux over the interval
# `[times[i], times[i] + dt)`. When a boundary falls inside that interval,
# the sample value is the area-weighted mix of the two phases' fluxes.
# This makes the rectangle-rule integral `sum(flux) × dt` equal the total
# injected iodine exactly, independent of whether each phase's duration is
# a multiple of dt. Physically the injector also doesn't instantaneously
# switch rates — there's always a finite ramp on the order of one
# injector control cycle, which this dt-wide blending coarsely models.

function _rectangular_phases_profile(phases, dt::Float64, t_max::Float64)
    times = collect(0.0:dt:t_max)
    flux  = zeros(length(times))
    cum_t = 0.0
    @inbounds for ph in phases
        t_lo = cum_t
        t_hi = cum_t + ph.volume_ml / ph.rate_ml_s
        cum_t = t_hi
        phase_flux = ph.rate_ml_s * ph.concentration_mgI_ml
        (t_hi <= t_lo || phase_flux <= 0.0) && continue
        for i in eachindex(times)
            t = times[i]
            overlap_lo = max(t, t_lo)
            overlap_hi = min(t + dt, t_hi)
            if overlap_hi > overlap_lo
                flux[i] += phase_flux * (overlap_hi - overlap_lo) / dt
            end
        end
    end
    return times, flux
end

# ── injection_profile dispatch ─────────────────────────────────────────────

function injection_profile(p::UniphaseNoChaser; dt::Float64=0.1, t_max::Float64=60.0)
    vol = p.weight_kg * p.contrast_volume_per_kg
    phases = ((volume_ml=vol, rate_ml_s=p.injection_rate_ml_s,
               concentration_mgI_ml=p.contrast_concentration_mgI_ml),)
    return _rectangular_phases_profile(phases, dt, t_max)
end

function injection_profile(p::UniphaseWithChaser; dt::Float64=0.1, t_max::Float64=60.0)
    contrast_vol = p.weight_kg * p.contrast_volume_per_kg
    chaser_vol   = p.weight_kg * p.chaser_volume_per_kg
    chaser_conc  = p.contrast_concentration_mgI_ml * p.chaser_dilution
    phases = (
        (volume_ml=contrast_vol, rate_ml_s=p.injection_rate_ml_s,
         concentration_mgI_ml=p.contrast_concentration_mgI_ml),
        (volume_ml=chaser_vol,   rate_ml_s=p.chaser_rate_ml_s,
         concentration_mgI_ml=chaser_conc),
    )
    return _rectangular_phases_profile(phases, dt, t_max)
end

function injection_profile(p::BiphaseNoChaser; dt::Float64=0.1, t_max::Float64=60.0)
    v1 = p.weight_kg * p.phase1_volume_per_kg
    v2 = p.weight_kg * p.phase2_volume_per_kg
    phases = (
        (volume_ml=v1, rate_ml_s=p.phase1_rate_ml_s,
         concentration_mgI_ml=p.contrast_concentration_mgI_ml),
        (volume_ml=v2, rate_ml_s=p.phase2_rate_ml_s,
         concentration_mgI_ml=p.contrast_concentration_mgI_ml),
    )
    return _rectangular_phases_profile(phases, dt, t_max)
end

function injection_profile(p::BiphaseWithChaser; dt::Float64=0.1, t_max::Float64=60.0)
    v1 = p.weight_kg * p.phase1_volume_per_kg
    v2 = p.weight_kg * p.phase2_volume_per_kg
    vc = p.weight_kg * p.chaser_volume_per_kg
    cc = p.contrast_concentration_mgI_ml * p.chaser_dilution
    phases = (
        (volume_ml=v1, rate_ml_s=p.phase1_rate_ml_s,
         concentration_mgI_ml=p.contrast_concentration_mgI_ml),
        (volume_ml=v2, rate_ml_s=p.phase2_rate_ml_s,
         concentration_mgI_ml=p.contrast_concentration_mgI_ml),
        (volume_ml=vc, rate_ml_s=p.chaser_rate_ml_s,
         concentration_mgI_ml=cc),
    )
    return _rectangular_phases_profile(phases, dt, t_max)
end

# ── Central-circulation gamma-variate kernel ──────────────────────────────

"""
    central_transit_kernel(times, peak_s, dispersion_s) -> Vector{Float64}

Right-skewed gamma-variate impulse response of the central circulation
(RV → lung → LV → aorta), normalized to ∫g dt = 1 (mass-preserving).

  g(t) = A × (t/T)^α × exp(α × (1 - t/T))     for t > 0
  where T = `peak_s`, α = (peak_s / dispersion_s)²

Larger α → narrower peak. Default `(12.0, 3.0)` gives α ≈ 16
(reasonably tight peak around 12 s). The gamma shape captures the long
recirculation tail seen in clinical AIF curves.
"""
function central_transit_kernel(times::Vector{Float64}, peak_s::Float64, dispersion_s::Float64)
    g = zeros(length(times))
    α = (peak_s / dispersion_s)^2
    @inbounds for (i, t) in enumerate(times)
        if t > 0
            tp = t / peak_s
            g[i] = tp^α * exp(α * (1.0 - tp))
        end
    end
    dt = length(times) > 1 ? times[2] - times[1] : 1.0
    Z = sum(g) * dt
    Z > 0 && (g ./= Z)
    return g
end

# ── Causal discrete convolution: (x ⊛ k)[i] = Σ_{j≤i} x[i-j+1] k[j] dt ─────

function _convolve_causal(x::Vector{Float64}, k::Vector{Float64}, dt::Float64)
    n = length(x)
    nk = length(k)
    y = zeros(n)
    @inbounds for i in 1:n
        s = 0.0
        j_max = min(nk, i)
        for j in 1:j_max
            s += x[i - j + 1] * k[j]
        end
        y[i] = s * dt
    end
    return y
end

# ── Top-level: protocol + physiology → AIF at aorta root ──────────────────

"""
    protocol_to_aif(protocol, physiology; dt=0.1, t_max=60.0)
        -> (times::Vector{Float64}, aif_mgI_per_ml::Vector{Float64})

Synthesize the AIF (arterial input function) at the aorta root.

Steps:
1. `injection_profile(protocol)` → mgI/s into peripheral vein
2. Convolve with `central_transit_kernel(physiology)` → mgI/s at aorta
3. ÷ `cardiac_output_ml_s` → mgI/mL at aorta

Mass conservation: ∫ aif dt = (total injected iodine) / cardiac_output.
"""
function protocol_to_aif(protocol::AbstractInjectionProtocol,
                          physiology::PatientPhysiology;
                          dt::Float64=0.1, t_max::Float64=60.0)
    times, inj_flux = injection_profile(protocol; dt=dt, t_max=t_max)
    kernel = central_transit_kernel(times,
        physiology.central_transit_delay_s,
        physiology.central_transit_dispersion_s)
    arterial_flux  = _convolve_causal(inj_flux, kernel, dt)         # mgI/s at aorta
    aif_mgI_per_ml = arterial_flux ./ physiology.cardiac_output_ml_s
    return times, aif_mgI_per_ml
end
