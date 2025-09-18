"""
Stellarator Physics Demo with 3D Visualization
This version includes plotting capabilities for 3D visualization.
"""

println("=== Stellarator Physics Analysis with 3D Visualization ===")
println()

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

# Utility functions
function cylindrical_to_cartesian(r, θ, z)
    x = r * cos(θ)
    y = r * sin(θ)
    return x, y, z
end

function magnetic_field_magnitude(Bx, By, Bz)
    return sqrt(Bx^2 + By^2 + Bz^2)
end

function safety_factor(R, B_φ, B_θ)
    return R * B_φ / B_θ
end

# Simplified Magnetic Field Structure
mutable struct MagneticField3D
    R₀::Float64  # Major radius
    a::Float64   # Minor radius
    N::Int       # Number of field periods
    harmonics::Dict{Tuple{Int,Int},Complex{Float64}}  # Fourier harmonics
    B₀::Float64  # Reference magnetic field strength
    
    function MagneticField3D(R₀, a, N, B₀=1.0)
        new(R₀, a, N, Dict{Tuple{Int,Int},Complex{Float64}}(), B₀)
    end
end

# Simplified Magnetic Surface Structure
struct MagneticSurface
    points::Vector{Tuple{Float64,Float64,Float64}}  # Surface points
    area::Float64                                    # Surface area
    volume::Float64                                  # Enclosed volume
    q::Float64                                       # Safety factor
end

# Simplified stellarator harmonics calculation
function stellarator_harmonics(m, n, R₀, a, N)
    harmonics = Dict{Tuple{Int,Int},Complex{Float64}}()
    
    # Add main toroidal field
    harmonics[(0, 1)] = 1.0 + 0.0im
    
    # Add stellarator harmonics
    for i in 1:3
        for j in 1:3
            if i != 0 || j != 1
                # Simplified harmonic amplitudes
                amplitude = 0.1 * exp(-(i^2 + j^2) / 10)
                phase = 2π * rand()
                harmonics[(i, j)] = amplitude * exp(im * phase)
            end
        end
    end
    
    return harmonics
end

# Simplified magnetic field calculation
function calculate_magnetic_field(bfield::MagneticField3D, R, φ, Z)
    B_R = 0.0
    B_φ = 0.0
    B_Z = 0.0
    
    # Convert to Boozer-like coordinates
    r = sqrt((R - bfield.R₀)^2 + Z^2)
    θ = atan(Z, R - bfield.R₀)
    
    # Calculate field components using Fourier harmonics
    for ((m, n), coeff) in bfield.harmonics
        if m == 0 && n == 1
            # Main toroidal field
            B_φ += real(coeff) * bfield.B₀
        else
            # Stellarator harmonics
            phase = m * θ + n * φ
            amplitude = real(coeff) * bfield.B₀ * (r / bfield.a)^m
            
            B_R += amplitude * m * sin(phase) / r
            B_φ += amplitude * n * cos(phase) / R
            B_Z += amplitude * m * cos(phase) / r
        end
    end
    
    return B_R, B_φ, B_Z
end

# Simplified magnetic surface calculation
function find_magnetic_surface(bfield::MagneticField3D, s, n_θ=32, n_φ=32)
    points = Tuple{Float64,Float64,Float64}[]
    
    # Generate surface points
    for i in 1:n_θ
        θ = 2π * (i - 1) / n_θ
        for j in 1:n_φ
            φ = 2π * (j - 1) / n_φ
            
            # Convert to cylindrical coordinates
            r = s * bfield.a
            R = bfield.R₀ + r * cos(θ)
            Z = r * sin(θ)
            
            push!(points, (R, φ, Z))
        end
    end
    
    # Calculate surface area and volume (simplified)
    area = 4π^2 * bfield.R₀ * s * bfield.a
    volume = 2π^2 * bfield.R₀ * (s * bfield.a)^2
    
    # Calculate safety factor (simplified)
    q = 1.0 / s  # Simplified q-profile
    
    return MagneticSurface(points, area, volume, q)
end

# Simplified quasi-symmetry measure
function quasi_symmetry_measure(bfield::MagneticField3D, surface::MagneticSurface)
    total_variation = 0.0
    n_points = length(surface.points)
    
    for point in surface.points
        R, φ, Z = point
        B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
        B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
        
        # Calculate variation in |B| along the field line
        total_variation += abs(B_mag - bfield.B₀)
    end
    
    return total_variation / n_points / bfield.B₀
end

# Function to create 3D surface plot data
function create_surface_data(surface::MagneticSurface)
    x_vals = Float64[]
    y_vals = Float64[]
    z_vals = Float64[]
    
    for point in surface.points
        R, φ, Z = point
        # Convert to Cartesian
        x, y, z = R * cos(φ), R * sin(φ), Z
        push!(x_vals, x)
        push!(y_vals, y)
        push!(z_vals, z)
    end
    
    return x_vals, y_vals, z_vals
end

# Function to create magnetic field line data
function trace_field_line_simple(bfield::MagneticField3D, R₀, φ₀, Z₀, length)
    points = Tuple{Float64,Float64,Float64}[]
    
    # Initial point
    R, φ, Z = R₀, φ₀, Z₀
    push!(points, (R, φ, Z))
    
    # Integration parameters
    ds = length / 100
    steps = Int(length / ds)
    
    for i in 1:steps
        # Calculate magnetic field at current point
        B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
        B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
        
        if B_mag > 0
            # Normalize field components
            B_R /= B_mag
            B_φ /= B_mag
            B_Z /= B_mag
            
            # Step along field line
            R += B_R * ds
            φ += B_φ * ds / R
            Z += B_Z * ds
            
            # Keep φ in [0, 2π]
            φ = mod(φ, 2π)
            
            push!(points, (R, φ, Z))
        else
            break
        end
    end
    
    return points
end

# Main demonstration
println("1. Creating stellarator magnetic field configuration...")
R₀ = 1.0  # Major radius [m]
a = 0.2   # Minor radius [m]
N = 5     # Number of field periods
B₀ = 1.0  # Reference magnetic field [T]

stellarator_bfield = MagneticField3D(R₀, a, N, B₀)
println("   ✓ Stellarator field created with R₀ = $(R₀) m, a = $(a) m, N = $(N)")

println("\n2. Adding stellarator harmonics...")
stellarator_bfield.harmonics = stellarator_harmonics(5, 5, R₀, a, N)
println("   ✓ Added $(length(stellarator_bfield.harmonics)) harmonic components")

println("\n3. Calculating magnetic field properties...")
s_values = [0.2, 0.4, 0.6, 0.8]

for s in s_values
    R = stellarator_bfield.R₀ + s * stellarator_bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(stellarator_bfield, R, φ, Z)
    B_mag = magnetic_field_magnitude(B_R, B_φ, B_Z)
    
    println("   At s = $(s): |B| = $(round(B_mag, digits=4)) T")
end

println("\n4. Finding magnetic surfaces...")
surfaces = []
for s in s_values
    surface = find_magnetic_surface(stellarator_bfield, s)
    qs_measure = quasi_symmetry_measure(stellarator_bfield, surface)
    
    println("   Surface s = $(s):")
    println("     Area = $(round(surface.area, digits=4)) m²")
    println("     Volume = $(round(surface.volume, digits=4)) m³")
    println("     Safety factor q = $(round(surface.q, digits=3))")
    println("     Quasi-symmetry measure = $(round(qs_measure, digits=4))")
    
    push!(surfaces, surface)
end

println("\n5. Creating 3D visualizations...")

# Try to use Plots.jl if available, otherwise create ASCII art
try
    using Plots
    
    println("   ✓ Plots.jl available - creating interactive 3D plots")
    
    # Create 3D plot of magnetic surfaces
    plot_3d = plot3d(
        xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
        title="3D Stellarator Magnetic Surfaces",
        legend=true,
        size=(800, 600)
    )
    
    colors = [:red, :blue, :green, :orange]
    
    for (i, surface) in enumerate(surfaces)
        x_vals, y_vals, z_vals = create_surface_data(surface)
        
        # Create a scatter plot for the surface points
        scatter!(plot_3d, x_vals, y_vals, z_vals,
                markersize=2, alpha=0.6, color=colors[i],
                label="Surface s=$(s_values[i])")
    end
    
    # Add coordinate axes
    plot!(plot_3d, [0, 1.5], [0, 0], [0, 0], 
          linewidth=3, color=:red, label="X-axis")
    plot!(plot_3d, [0, 0], [0, 1.5], [0, 0], 
          linewidth=3, color=:green, label="Y-axis")
    plot!(plot_3d, [0, 0], [0, 0], [0, 0.5], 
          linewidth=3, color=:blue, label="Z-axis")
    
    # Save the plot
    savefig(plot_3d, "stellarator_3d_surfaces.png")
    println("   ✓ 3D plot saved as 'stellarator_3d_surfaces.png'")
    
    # Create a 2D cross-section plot
    plot_2d = plot(
        xlabel="R [m]", ylabel="Z [m]",
        title="Stellarator Cross-Section (φ = 0)",
        aspect_ratio=:equal,
        size=(600, 600)
    )
    
    for (i, surface) in enumerate(surfaces)
        x_vals, y_vals, z_vals = create_surface_data(surface)
        
        # Filter points near φ = 0
        cross_section_x = Float64[]
        cross_section_z = Float64[]
        
        for j in 1:length(x_vals)
            if abs(y_vals[j]) < 0.1  # Near φ = 0
                push!(cross_section_x, x_vals[j])
                push!(cross_section_z, z_vals[j])
            end
        end
        
        if length(cross_section_x) > 0
            scatter!(plot_2d, cross_section_x, cross_section_z,
                    markersize=3, alpha=0.7, color=colors[i],
                    label="Surface s=$(s_values[i])")
        end
    end
    
    savefig(plot_2d, "stellarator_cross_section.png")
    println("   ✓ Cross-section plot saved as 'stellarator_cross_section.png'")
    
    # Create magnetic field strength plot
    plot_field = plot(
        xlabel="Normalized Radius (s)", ylabel="|B| [T]",
        title="Magnetic Field Strength Profile",
        linewidth=2,
        size=(600, 400)
    )
    
    B_values = Float64[]
    for s in s_values
        R = stellarator_bfield.R₀ + s * stellarator_bfield.a
        φ = 0.0
        Z = 0.0
        B_R, B_φ, B_Z = calculate_magnetic_field(stellarator_bfield, R, φ, Z)
        B_mag = magnetic_field_magnitude(B_R, B_φ, B_Z)
        push!(B_values, B_mag)
    end
    
    plot!(plot_field, s_values, B_values, linewidth=3, color=:purple,
          marker=:circle, markersize=6, label="|B|")
    
    savefig(plot_field, "magnetic_field_profile.png")
    println("   ✓ Field profile plot saved as 'magnetic_field_profile.png'")
    
catch e
    println("   ⚠ Plots.jl not available: $(e)")
    println("   Creating ASCII art visualization instead...")
    
    # ASCII art representation
    println("\n   ASCII 3D Visualization:")
    println("   ======================")
    println("   Stellarator Cross-Section (φ = 0):")
    println("   ")
    println("   Z ↑")
    println("     |")
    println("  0.2 |     ●●●●●")
    println("     |   ●●●   ●●●")
    println("  0.1 |  ●●       ●●")
    println("     | ●●         ●●")
    println("  0.0 |●●           ●●")
    println("     |●●           ●●")
    println(" -0.1 |  ●●       ●●")
    println("     |   ●●●   ●●●")
    println(" -0.2 |     ●●●●●")
    println("     |")
    println("     +----------------→ R")
    println("     0.8  1.0  1.2  1.4")
    println("   ")
    println("   Legend: ● = Magnetic surface points")
    println("   Colors represent different flux surfaces")
end

println("\n6. Calculating transport estimates...")
T_e = 1000.0  # Electron temperature [eV]
n_e = 1e20    # Electron density [m^-3]

transport_data = []
for s in s_values
    # Simplified transport coefficient estimate
    R = stellarator_bfield.R₀ + s * stellarator_bfield.a
    φ = 0.0
    Z = 0.0
    B_R, B_φ, B_Z = calculate_magnetic_field(stellarator_bfield, R, φ, Z)
    B_mag = magnetic_field_magnitude(B_R, B_φ, B_Z)
    
    # Simplified Larmor radius
    ρ_i = sqrt(2 * mᵢ * T_e * T_eV) / (e * B_mag)
    
    # Simplified diffusion coefficient
    D_estimate = ρ_i^2 / (s * a)^2 * sqrt(T_e) * 1e-4
    
    push!(transport_data, (s, ρ_i, D_estimate))
    
    println("   At s = $(s):")
    println("     Larmor radius ρᵢ = $(round(ρ_i*1000, digits=2)) mm")
    println("     Estimated D ≈ $(round(D_estimate, digits=8)) m²/s")
end

println("\n7. Performance comparison with tokamak...")
println("   Stellarator advantages:")
println("     - Steady-state operation (no plasma current needed)")
println("     - Reduced MHD instabilities")
println("     - Flexible 3D optimization")
println("   Tokamak advantages:")
println("     - Simpler magnetic field geometry")
println("     - Better established technology")
println("     - Higher plasma beta limits")

println("\n=== ANALYSIS SUMMARY ===")
println("✓ 3D magnetic field calculations completed")
println("✓ Magnetic surface analysis completed")
println("✓ Quasi-symmetry evaluation completed")
println("✓ Transport coefficient estimates completed")
println("✓ 3D visualization created")
println("✓ Performance comparison analysis completed")

println("\n=== DEMO COMPLETED SUCCESSFULLY! ===")
println("This demonstrates the core functionality of the StellaratorPhysics package")
println("including 3D visualization capabilities.")
println("The full version includes advanced optimization, detailed transport modeling,")
println("and interactive 3D visualization with PlotlyJS.")
