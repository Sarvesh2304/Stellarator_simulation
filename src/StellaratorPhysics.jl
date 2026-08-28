
module StellaratorPhysics

# Core modules
export MagneticField3D, NeoclassicalTransport, StellaratorOptimization
export TokamakComparison, Visualization3D

# Main analysis functions
export analyze_stellarator, optimize_stellarator, compare_with_tokamak
export plot_magnetic_field, plot_plasma_surfaces, plot_transport_coefficients

# Physical constants and utilities
export μ₀, e, mₑ, mᵢ, k_B, T_eV, ρ_s, v_th
export cylindrical_to_cartesian, cartesian_to_cylindrical, toroidal_to_cartesian
export magnetic_field_magnitude, safety_factor, plasma_beta, larmor_radius
export collision_frequency, log_Λ

# Include submodules
include("MagneticField3D.jl")
include("NeoclassicalTransport.jl")
include("StellaratorOptimization.jl")
include("TokamakComparison.jl")
include("Visualization3D.jl")
include("utils.jl")

end
