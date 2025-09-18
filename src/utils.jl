"""
Physical constants and utility functions for stellarator physics analysis
"""

# Physical constants
const μ₀ = 4π * 1e-7  # Permeability of free space [H/m]
const e = 1.602176634e-19  # Elementary charge [C]
const mₑ = 9.1093837015e-31  # Electron mass [kg]
const mᵢ = 1.67262192369e-27  # Proton mass [kg]
const k_B = 1.380649e-23  # Boltzmann constant [J/K]
const T_eV = e / k_B  # Temperature conversion factor [K/eV]

# Derived constants
const ρ_s = sqrt(mᵢ * T_eV / (e^2 * μ₀))  # Sound gyroradius scale
const v_th = sqrt(2 * e * T_eV / mₑ)  # Thermal velocity scale

"""
    cylindrical_to_cartesian(r, θ, z)

Convert cylindrical coordinates to Cartesian coordinates.
"""
function cylindrical_to_cartesian(r, θ, z)
    x = r * cos(θ)
    y = r * sin(θ)
    return x, y, z
end

"""
    cartesian_to_cylindrical(x, y, z)

Convert Cartesian coordinates to cylindrical coordinates.
"""
function cartesian_to_cylindrical(x, y, z)
    r = sqrt(x^2 + y^2)
    θ = atan(y, x)
    return r, θ, z
end

"""
    toroidal_to_cartesian(R, φ, Z)

Convert toroidal coordinates to Cartesian coordinates.
"""
function toroidal_to_cartesian(R, φ, Z)
    x = R * cos(φ)
    y = R * sin(φ)
    z = Z
    return x, y, z
end

"""
    magnetic_field_magnitude(Bx, By, Bz)

Calculate the magnitude of a 3D magnetic field vector.
"""
function magnetic_field_magnitude(Bx, By, Bz)
    return sqrt(Bx^2 + By^2 + Bz^2)
end

"""
    safety_factor(R, B_φ, B_θ)

Calculate the safety factor q for a given magnetic field configuration.
"""
function safety_factor(R, B_φ, B_θ)
    return R * B_φ / B_θ
end

"""
    plasma_beta(pressure, B_mag)

Calculate the plasma beta parameter.
"""
function plasma_beta(pressure, B_mag)
    return 2 * μ₀ * pressure / B_mag^2
end

"""
    larmor_radius(mass, temperature, B_mag)

Calculate the Larmor radius for a particle in a magnetic field.
"""
function larmor_radius(mass, temperature, B_mag)
    return sqrt(2 * mass * temperature) / (e * B_mag)
end

"""
    collision_frequency(n, T, Z)

Calculate the collision frequency for plasma particles.
"""
function collision_frequency(n, T, Z)
    # Simplified collision frequency formula
    return 1.6e-19 * n * Z^2 * e^4 * log(Λ) / (4π * ε₀^2 * mₑ^2 * v_th^3)
end

"""
    log_Λ(n, T)

Calculate the Coulomb logarithm.
"""
function log_Λ(n, T)
    return 23 - log(sqrt(n) * T^(-3/2))
end
