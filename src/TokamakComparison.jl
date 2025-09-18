"""
Tokamak Comparison Module

This module provides comparison capabilities between stellarator and tokamak
configurations, including:
- Performance metrics comparison
- Transport coefficient comparison
- Stability analysis comparison
- Confinement time scaling
"""

using LinearAlgebra
using Statistics

export TokamakComparison, ComparisonResult
export compare_performance, compare_transport, compare_stability
export confinement_scaling, beta_limits

"""
    struct TokamakComparison

Main structure for tokamak comparison calculations.
"""
mutable struct TokamakComparison
    stellarator_bfield::MagneticField3D
    tokamak_bfield::MagneticField3D
    stellarator_transport::NeoclassicalTransport
    tokamak_transport::NeoclassicalTransport
    
    function TokamakComparison(stellarator_bfield, tokamak_bfield, 
                              stellarator_transport, tokamak_transport)
        new(stellarator_bfield, tokamak_bfield, stellarator_transport, tokamak_transport)
    end
end

"""
    struct ComparisonResult

Stores the results of a stellarator-tokamak comparison.
"""
struct ComparisonResult
    metric::String
    stellarator_value::Float64
    tokamak_value::Float64
    ratio::Float64
    advantage::String
end

"""
    create_tokamak_field(R₀, a, q₀, q_a, B₀)

Create a simplified tokamak magnetic field configuration.
"""
function create_tokamak_field(R₀, a, q₀, q_a, B₀)
    # Create tokamak field with simplified q-profile
    bfield = MagneticField3D(R₀, a, 1, B₀)  # N=1 for tokamak
    
    # Add tokamak-specific harmonics
    harmonics = Dict{Tuple{Int,Int},Complex{Float64}}()
    
    # Main toroidal field
    harmonics[(0, 1)] = 1.0 + 0.0im
    
    # Poloidal field (simplified)
    harmonics[(1, 0)] = 0.1 + 0.0im
    
    bfield.harmonics = harmonics
    return bfield
end

"""
    compare_performance(comparison::TokamakComparison, s)

Compare performance metrics between stellarator and tokamak.
"""
function compare_performance(comparison::TokamakComparison, s)
    results = ComparisonResult[]
    
    # 1. Safety factor comparison
    stellarator_surface = find_magnetic_surface(comparison.stellarator_bfield, s)
    tokamak_surface = find_magnetic_surface(comparison.tokamak_bfield, s)
    
    q_stellarator = stellarator_surface.q
    q_tokamak = tokamak_surface.q
    
    push!(results, ComparisonResult(
        "Safety Factor",
        q_stellarator,
        q_tokamak,
        q_stellarator / q_tokamak,
        q_stellarator > q_tokamak ? "Stellarator" : "Tokamak"
    ))
    
    # 2. Magnetic field strength comparison
    R = comparison.stellarator_bfield.R₀ + s * comparison.stellarator_bfield.a
    φ = 0.0
    Z = 0.0
    
    B_R_s, B_φ_s, B_Z_s = calculate_magnetic_field(comparison.stellarator_bfield, R, φ, Z)
    B_mag_s = sqrt(B_R_s^2 + B_φ_s^2 + B_Z_s^2)
    
    B_R_t, B_φ_t, B_Z_t = calculate_magnetic_field(comparison.tokamak_bfield, R, φ, Z)
    B_mag_t = sqrt(B_R_t^2 + B_φ_t^2 + B_Z_t^2)
    
    push!(results, ComparisonResult(
        "Magnetic Field Strength",
        B_mag_s,
        B_mag_t,
        B_mag_s / B_mag_t,
        B_mag_s > B_mag_t ? "Stellarator" : "Tokamak"
    ))
    
    # 3. Plasma volume comparison
    push!(results, ComparisonResult(
        "Plasma Volume",
        stellarator_surface.volume,
        tokamak_surface.volume,
        stellarator_surface.volume / tokamak_surface.volume,
        stellarator_surface.volume > tokamak_surface.volume ? "Stellarator" : "Tokamak"
    ))
    
    return results
end

"""
    compare_transport(comparison::TokamakComparison, s)

Compare transport coefficients between stellarator and tokamak.
"""
function compare_transport(comparison::TokamakComparison, s)
    results = ComparisonResult[]
    
    # Get transport coefficients
    stellarator_coeffs = calculate_transport_coefficients(comparison.stellarator_transport, s)
    tokamak_coeffs = calculate_transport_coefficients(comparison.tokamak_transport, s)
    
    # 1. Particle diffusion coefficient
    push!(results, ComparisonResult(
        "Particle Diffusion (D_11)",
        stellarator_coeffs.D_11,
        tokamak_coeffs.D_11,
        stellarator_coeffs.D_11 / tokamak_coeffs.D_11,
        stellarator_coeffs.D_11 < tokamak_coeffs.D_11 ? "Stellarator" : "Tokamak"
    ))
    
    # 2. Heat diffusion coefficient
    push!(results, ComparisonResult(
        "Heat Diffusion (D_22)",
        stellarator_coeffs.D_22,
        tokamak_coeffs.D_22,
        stellarator_coeffs.D_22 / tokamak_coeffs.D_22,
        stellarator_coeffs.D_22 < tokamak_coeffs.D_22 ? "Stellarator" : "Tokamak"
    ))
    
    # 3. Effective thermal diffusivity
    push!(results, ComparisonResult(
        "Thermal Diffusivity (χ_eff)",
        stellarator_coeffs.χ_eff,
        tokamak_coeffs.χ_eff,
        stellarator_coeffs.χ_eff / tokamak_coeffs.χ_eff,
        stellarator_coeffs.χ_eff < tokamak_coeffs.χ_eff ? "Stellarator" : "Tokamak"
    ))
    
    # 4. Bootstrap current
    j_bs_s = bootstrap_current(comparison.stellarator_transport, s)
    j_bs_t = bootstrap_current(comparison.tokamak_transport, s)
    
    push!(results, ComparisonResult(
        "Bootstrap Current",
        j_bs_s,
        j_bs_t,
        j_bs_s / j_bs_t,
        j_bs_s > j_bs_t ? "Stellarator" : "Tokamak"
    ))
    
    return results
end

"""
    compare_stability(comparison::TokamakComparison, s)

Compare stability properties between stellarator and tokamak.
"""
function compare_stability(comparison::TokamakComparison, s)
    results = ComparisonResult[]
    
    # 1. Magnetic well depth
    well_depth_s = magnetic_well_depth(comparison.stellarator_bfield, 0.1, 0.9)
    well_depth_t = magnetic_well_depth(comparison.tokamak_bfield, 0.1, 0.9)
    
    push!(results, ComparisonResult(
        "Magnetic Well Depth",
        well_depth_s,
        well_depth_t,
        well_depth_s / well_depth_t,
        well_depth_s > well_depth_t ? "Stellarator" : "Tokamak"
    ))
    
    # 2. Beta limit (simplified)
    beta_limit_s = calculate_beta_limit(comparison.stellarator_bfield, s)
    beta_limit_t = calculate_beta_limit(comparison.tokamak_bfield, s)
    
    push!(results, ComparisonResult(
        "Beta Limit",
        beta_limit_s,
        beta_limit_t,
        beta_limit_s / beta_limit_t,
        beta_limit_s > beta_limit_t ? "Stellarator" : "Tokamak"
    ))
    
    # 3. Mercier stability (simplified)
    mercier_s = calculate_mercier_criterion(comparison.stellarator_bfield, s)
    mercier_t = calculate_mercier_criterion(comparison.tokamak_bfield, s)
    
    push!(results, ComparisonResult(
        "Mercier Criterion",
        mercier_s,
        mercier_t,
        mercier_s / mercier_t,
        mercier_s > mercier_t ? "Stellarator" : "Tokamak"
    ))
    
    return results
end

"""
    calculate_beta_limit(bfield::MagneticField3D, s)

Calculate the beta limit for a given magnetic field configuration.
"""
function calculate_beta_limit(bfield::MagneticField3D, s)
    # Simplified beta limit calculation
    # Real implementations would use more sophisticated stability analysis
    
    R = bfield.R₀ + s * bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
    B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
    
    # Simplified beta limit based on magnetic field strength
    beta_limit = 0.1 * (B_mag / bfield.B₀)^2
    
    return beta_limit
end

"""
    calculate_mercier_criterion(bfield::MagneticField3D, s)

Calculate the Mercier stability criterion.
"""
function calculate_mercier_criterion(bfield::MagneticField3D, s)
    # Simplified Mercier criterion calculation
    # Real implementations would use more sophisticated analysis
    
    R = bfield.R₀ + s * bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
    B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
    
    # Simplified Mercier criterion
    mercier = 1.0 - 0.5 * s^2  # Simplified model
    
    return mercier
end

"""
    confinement_scaling(comparison::TokamakComparison, s)

Compare confinement time scaling between stellarator and tokamak.
"""
function confinement_scaling(comparison::TokamakComparison, s)
    results = ComparisonResult[]
    
    # Calculate confinement times
    τ_E_s = neoclassical_energy_confinement(comparison.stellarator_transport, s)
    τ_E_t = neoclassical_energy_confinement(comparison.tokamak_transport, s)
    
    push!(results, ComparisonResult(
        "Energy Confinement Time",
        τ_E_s,
        τ_E_t,
        τ_E_s / τ_E_t,
        τ_E_s > τ_E_t ? "Stellarator" : "Tokamak"
    ))
    
    # Calculate scaling with plasma parameters
    # This is a simplified scaling analysis
    scaling_s = calculate_scaling_law(comparison.stellarator_bfield, s)
    scaling_t = calculate_scaling_law(comparison.tokamak_bfield, s)
    
    push!(results, ComparisonResult(
        "Scaling Law Exponent",
        scaling_s,
        scaling_t,
        scaling_s / scaling_t,
        scaling_s > scaling_t ? "Stellarator" : "Tokamak"
    ))
    
    return results
end

"""
    calculate_scaling_law(bfield::MagneticField3D, s)

Calculate the scaling law exponent for confinement time.
"""
function calculate_scaling_law(bfield::MagneticField3D, s)
    # Simplified scaling law calculation
    # Real implementations would use more sophisticated analysis
    
    R = bfield.R₀ + s * bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
    B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
    
    # Simplified scaling exponent
    scaling_exponent = 1.0 + 0.1 * s  # Simplified model
    
    return scaling_exponent
end

"""
    beta_limits(comparison::TokamakComparison, s)

Compare beta limits between stellarator and tokamak.
"""
function beta_limits(comparison::TokamakComparison, s)
    results = ComparisonResult[]
    
    # Calculate beta limits
    beta_limit_s = calculate_beta_limit(comparison.stellarator_bfield, s)
    beta_limit_t = calculate_beta_limit(comparison.tokamak_bfield, s)
    
    push!(results, ComparisonResult(
        "Beta Limit",
        beta_limit_s,
        beta_limit_t,
        beta_limit_s / beta_limit_t,
        beta_limit_s > beta_limit_t ? "Stellarator" : "Tokamak"
    ))
    
    # Calculate beta limits for different modes
    beta_limit_s_ballooning = calculate_ballooning_beta_limit(comparison.stellarator_bfield, s)
    beta_limit_t_ballooning = calculate_ballooning_beta_limit(comparison.tokamak_bfield, s)
    
    push!(results, ComparisonResult(
        "Ballooning Beta Limit",
        beta_limit_s_ballooning,
        beta_limit_t_ballooning,
        beta_limit_s_ballooning / beta_limit_t_ballooning,
        beta_limit_s_ballooning > beta_limit_t_ballooning ? "Stellarator" : "Tokamak"
    ))
    
    return results
end

"""
    calculate_ballooning_beta_limit(bfield::MagneticField3D, s)

Calculate the ballooning beta limit.
"""
function calculate_ballooning_beta_limit(bfield::MagneticField3D, s)
    # Simplified ballooning beta limit calculation
    # Real implementations would use more sophisticated analysis
    
    R = bfield.R₀ + s * bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
    B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
    
    # Simplified ballooning beta limit
    beta_limit = 0.05 * (B_mag / bfield.B₀)^2
    
    return beta_limit
end

"""
    comprehensive_comparison(comparison::TokamakComparison, s)

Perform a comprehensive comparison between stellarator and tokamak.
"""
function comprehensive_comparison(comparison::TokamakComparison, s)
    all_results = Dict{String,Vector{ComparisonResult}}()
    
    all_results["Performance"] = compare_performance(comparison, s)
    all_results["Transport"] = compare_transport(comparison, s)
    all_results["Stability"] = compare_stability(comparison, s)
    all_results["Confinement"] = confinement_scaling(comparison, s)
    all_results["Beta Limits"] = beta_limits(comparison, s)
    
    return all_results
end
