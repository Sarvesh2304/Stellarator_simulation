"""
Neoclassical Transport in Stellarators

This module implements neoclassical transport calculations specific to stellarators,
including:
- Monoenergetic transport coefficients
- Bootstrap current calculations
- Radial electric field effects
- Stellarator-specific transport regimes
"""

using LinearAlgebra
using SpecialFunctions

export NeoclassicalTransport, TransportCoefficients
export calculate_transport_coefficients, bootstrap_current
export radial_electric_field, neoclassical_energy_confinement

"""
    struct TransportCoefficients

Stores neoclassical transport coefficients for stellarators.
"""
struct TransportCoefficients
    D_11::Float64  # Particle diffusion coefficient
    D_12::Float64  # Cross-diffusion coefficient
    D_21::Float64  # Cross-diffusion coefficient
    D_22::Float64  # Heat diffusion coefficient
    χ_eff::Float64 # Effective thermal diffusivity
    ν_eff::Float64 # Effective viscosity
end

"""
    struct NeoclassicalTransport

Main structure for neoclassical transport calculations.
"""
mutable struct NeoclassicalTransport
    bfield::MagneticField3D
    T_e::Float64      # Electron temperature [eV]
    T_i::Float64      # Ion temperature [eV]
    n_e::Float64      # Electron density [m^-3]
    n_i::Float64      # Ion density [m^-3]
    Z_eff::Float64    # Effective charge
    ν_ei::Float64     # Electron-ion collision frequency
    ν_ii::Float64     # Ion-ion collision frequency
    
    function NeoclassicalTransport(bfield::MagneticField3D, T_e, T_i, n_e, n_i, Z_eff=1.0)
        # Calculate collision frequencies
        ν_ei = collision_frequency(n_e, T_e, Z_eff)
        ν_ii = collision_frequency(n_i, T_i, 1.0)
        
        new(bfield, T_e, T_i, n_e, n_i, Z_eff, ν_ei, ν_ii)
    end
end

"""
    collision_frequency(n, T, Z)

Calculate collision frequency for plasma particles.
"""
function collision_frequency(n, T, Z)
    # Coulomb collision frequency
    Λ = log_Λ(n, T)
    return 4π * n * Z^2 * e^4 * Λ / (3 * mₑ^2 * v_th^3 * (4π * ε₀)^2)
end

"""
    log_Λ(n, T)

Calculate the Coulomb logarithm.
"""
function log_Λ(n, T)
    return 23 - log(sqrt(n) * T^(-3/2))
end

"""
    monoenergetic_coefficients(transport::NeoclassicalTransport, s, energy)

Calculate monoenergetic transport coefficients for a given energy.
"""
function monoenergetic_coefficients(transport::NeoclassicalTransport, s, energy)
    # Get magnetic field properties at radius s
    R = transport.bfield.R₀ + s * transport.bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(transport.bfield, R, φ, Z)
    B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
    
    # Calculate Larmor radius
    ρ = larmor_radius(mᵢ, energy, B_mag)
    
    # Calculate effective collision frequency
    ν_eff = transport.ν_ii * (energy / transport.T_i)^(-3/2)
    
    # Calculate transport coefficients using stellarator theory
    # This is a simplified model - real implementations would use more sophisticated calculations
    
    # Particle diffusion coefficient
    D_11 = ρ^2 * ν_eff * (1 + 0.5 * s^2)
    
    # Cross-diffusion coefficient
    D_12 = -ρ^2 * ν_eff * s * 0.3
    
    # Heat diffusion coefficient
    D_22 = ρ^2 * ν_eff * (1 + 0.8 * s^2)
    
    # Effective thermal diffusivity
    χ_eff = D_22 / transport.n_i
    
    # Effective viscosity
    ν_eff_visc = ρ^2 * ν_eff * 0.1
    
    return TransportCoefficients(D_11, D_12, D_12, D_22, χ_eff, ν_eff_visc)
end

"""
    calculate_transport_coefficients(transport::NeoclassicalTransport, s)

Calculate integrated transport coefficients for a given radius.
"""
function calculate_transport_coefficients(transport::NeoclassicalTransport, s)
    # Energy integration for transport coefficients
    n_energy = 50
    energy_max = 5 * max(transport.T_e, transport.T_i)
    energies = range(0.1 * energy_max, energy_max, length=n_energy)
    
    # Weight function for energy integration
    weight_e = exp.(-energies / transport.T_e) .* sqrt(energies)
    weight_i = exp.(-energies / transport.T_i) .* sqrt(energies)
    
    # Initialize integrated coefficients
    D_11_int = 0.0
    D_12_int = 0.0
    D_22_int = 0.0
    χ_eff_int = 0.0
    ν_eff_int = 0.0
    
    # Integrate over energy
    for (i, energy) in enumerate(energies)
        coeffs = monoenergetic_coefficients(transport, s, energy)
        
        # Weight by Maxwellian distribution
        weight = weight_e[i] + weight_i[i]
        
        D_11_int += coeffs.D_11 * weight
        D_12_int += coeffs.D_12 * weight
        D_22_int += coeffs.D_22 * weight
        χ_eff_int += coeffs.χ_eff * weight
        ν_eff_int += coeffs.ν_eff * weight
    end
    
    # Normalize
    total_weight = sum(weight_e) + sum(weight_i)
    D_11_int /= total_weight
    D_12_int /= total_weight
    D_22_int /= total_weight
    χ_eff_int /= total_weight
    ν_eff_int /= total_weight
    
    return TransportCoefficients(D_11_int, D_12_int, D_12_int, D_22_int, χ_eff_int, ν_eff_int)
end

"""
    bootstrap_current(transport::NeoclassicalTransport, s)

Calculate the bootstrap current density at a given radius.
"""
function bootstrap_current(transport::NeoclassicalTransport, s)
    # Get transport coefficients
    coeffs = calculate_transport_coefficients(transport, s)
    
    # Calculate pressure gradient
    # This is a simplified model - real implementations would use actual pressure profiles
    p_grad = transport.n_e * transport.T_e / (transport.bfield.a * s)
    
    # Calculate bootstrap current density
    # Using simplified stellarator bootstrap current formula
    j_bs = -coeffs.D_12 * p_grad / (e * transport.bfield.B₀)
    
    return j_bs
end

"""
    radial_electric_field(transport::NeoclassicalTransport, s)

Calculate the radial electric field at a given radius.
"""
function radial_electric_field(transport::NeoclassicalTransport, s)
    # Get transport coefficients
    coeffs = calculate_transport_coefficients(transport, s)
    
    # Calculate density gradient
    # This is a simplified model
    n_grad = transport.n_e / (transport.bfield.a * s)
    
    # Calculate radial electric field
    # Using ambipolarity condition
    E_r = -coeffs.D_12 * n_grad / (e * transport.n_e)
    
    return E_r
end

"""
    neoclassical_energy_confinement(transport::NeoclassicalTransport, s)

Calculate the neoclassical energy confinement time.
"""
function neoclassical_energy_confinement(transport::NeoclassicalTransport, s)
    # Get transport coefficients
    coeffs = calculate_transport_coefficients(transport, s)
    
    # Calculate energy confinement time
    # Using simplified scaling
    τ_E = (transport.bfield.a * s)^2 / (2 * coeffs.χ_eff)
    
    return τ_E
end

"""
    stellarator_transport_regime(transport::NeoclassicalTransport, s)

Determine the dominant transport regime at a given radius.
"""
function stellarator_transport_regime(transport::NeoclassicalTransport, s)
    # Calculate relevant parameters
    R = transport.bfield.R₀ + s * transport.bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(transport.bfield, R, φ, Z)
    B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
    
    # Calculate Larmor radius
    ρ_i = larmor_radius(mᵢ, transport.T_i, B_mag)
    
    # Calculate collision frequency
    ν_ii = transport.ν_ii
    
    # Calculate transport regime
    if ν_ii * ρ_i / v_th < 0.1
        return "collisionless"
    elseif ν_ii * ρ_i / v_th < 1.0
        return "intermediate"
    else
        return "collisional"
    end
end

"""
    transport_scaling(transport::NeoclassicalTransport, s)

Calculate transport scaling with plasma parameters.
"""
function transport_scaling(transport::NeoclassicalTransport, s)
    # Get base transport coefficients
    coeffs = calculate_transport_coefficients(transport, s)
    
    # Calculate scaling with temperature
    T_scaling = (transport.T_e / 1000.0)^(3/2)
    
    # Calculate scaling with density
    n_scaling = (transport.n_e / 1e20)^(-1/2)
    
    # Calculate scaling with magnetic field
    B_scaling = (transport.bfield.B₀ / 1.0)^(-2)
    
    # Apply scaling
    D_11_scaled = coeffs.D_11 * T_scaling * n_scaling * B_scaling
    χ_eff_scaled = coeffs.χ_eff * T_scaling * n_scaling * B_scaling
    
    return TransportCoefficients(
        D_11_scaled, coeffs.D_12, coeffs.D_21, 
        coeffs.D_22 * T_scaling * n_scaling * B_scaling,
        χ_eff_scaled, coeffs.ν_eff
    )
end
