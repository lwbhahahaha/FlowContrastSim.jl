"""
    FlowContrastSim

A Julia package for hemodynamic flow simulation and dynamic contrast transport
in vascular trees.

# Quick Start
```julia
using FlowContrastSim

# From TOML config (full pipeline)
config = load_flow_config("configs/coronary_hyperemic.toml")
result = run_flow_simulation("path/to/tree_csvs", config; output_dir="output")

# Step-by-step
tree = load_tree("LAD", "lad_grown_segments.csv")
hemo = compute_hemodynamics(tree; target_flow_ml_min=242.0)
cr   = simulate_contrast(tree, hemo; dt=0.1, t_end=20.0)
```

Features:
- Poiseuille flow with Pries et al. (1992) non-Newtonian blood viscosity
- Terminal capillary bed resistance auto-calibration to match target flows
- Plug-flow contrast transport with dispersion broadening
- Interactive 3D contrast viewer (Plotly)

See also: [`compute_hemodynamics`](@ref), [`simulate_contrast`](@ref), [`run_flow_simulation`](@ref)
"""
module FlowContrastSim

using DelimitedFiles
using LinearAlgebra
using StaticArrays
using Statistics
using TOML

include("flow_tree.jl")
include("flow_config.jl")
include("hemodynamics.jl")
include("contrast_transport.jl")
include("contrast_viewer.jl")
include("run_simulation.jl")

# ── Types ──
export FlowTree,
       FlowConfig,
       HemodynamicsResult,
       ContrastResult

# ── Tree loading ──
export load_tree,
       load_trees

# ── Config ──
export load_flow_config

# ── Hemodynamics ──
export compute_hemodynamics,
       pries_viscosity_relative,
       apparent_viscosity

# ── Contrast simulation ──
export simulate_contrast,
       gamma_variate_input

# ── Visualization ──
export build_contrast_viewer

# ── Top-level pipeline ──
export run_flow_simulation

# ── Constants (for advanced users) ──
export BLOOD_VISCOSITY_PA_S,
       PLASMA_VISCOSITY_PA_S,
       DEFAULT_ROOT_PRESSURE_PA,
       DEFAULT_TERMINAL_PRESSURE_PA,
       DEFAULT_DISCHARGE_HEMATOCRIT

end # module FlowContrastSim
