"""
Basic test script for StellaratorPhysics package
"""

using Pkg
Pkg.activate("..")

using StellaratorPhysics
using Test

println("Testing StellaratorPhysics package...")

# Test 1: Create magnetic field
@testset "Magnetic Field Creation" begin
    bfield = MagneticField3D(1.0, 0.2, 5, 1.0)
    @test bfield.R₀ == 1.0
    @test bfield.a == 0.2
    @test bfield.N == 5
    @test bfield.B₀ == 1.0
    println("✓ Magnetic field creation test passed")
end

# Test 2: Calculate magnetic field
@testset "Magnetic Field Calculation" begin
    bfield = MagneticField3D(1.0, 0.2, 5, 1.0)
    bfield.harmonics = stellarator_harmonics(5, 5, 1.0, 0.2, 5)
    
    R, φ, Z = 1.2, 0.0, 0.1
    B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
    
    @test isfinite(B_R)
    @test isfinite(B_φ)
    @test isfinite(B_Z)
    println("✓ Magnetic field calculation test passed")
end

# Test 3: Create transport object
@testset "Transport Object Creation" begin
    bfield = MagneticField3D(1.0, 0.2, 5, 1.0)
    transport = NeoclassicalTransport(bfield, 1000.0, 1000.0, 1e20, 1e20, 1.0)
    
    @test transport.T_e == 1000.0
    @test transport.T_i == 1000.0
    @test transport.n_e == 1e20
    @test transport.n_i == 1e20
    println("✓ Transport object creation test passed")
end

# Test 4: Calculate transport coefficients
@testset "Transport Coefficient Calculation" begin
    bfield = MagneticField3D(1.0, 0.2, 5, 1.0)
    bfield.harmonics = stellarator_harmonics(5, 5, 1.0, 0.2, 5)
    transport = NeoclassicalTransport(bfield, 1000.0, 1000.0, 1e20, 1e20, 1.0)
    
    coeffs = calculate_transport_coefficients(transport, 0.5)
    
    @test isfinite(coeffs.D_11)
    @test isfinite(coeffs.D_22)
    @test isfinite(coeffs.χ_eff)
    println("✓ Transport coefficient calculation test passed")
end

# Test 5: Find magnetic surface
@testset "Magnetic Surface Calculation" begin
    bfield = MagneticField3D(1.0, 0.2, 5, 1.0)
    bfield.harmonics = stellarator_harmonics(5, 5, 1.0, 0.2, 5)
    
    surface = find_magnetic_surface(bfield, 0.5)
    
    @test length(surface.points) > 0
    @test isfinite(surface.area)
    @test isfinite(surface.volume)
    @test isfinite(surface.q)
    println("✓ Magnetic surface calculation test passed")
end

# Test 6: Utility functions
@testset "Utility Functions" begin
    # Test coordinate conversion
    x, y, z = cylindrical_to_cartesian(1.0, π/4, 0.0)
    @test isapprox(x, 1/sqrt(2), atol=1e-10)
    @test isapprox(y, 1/sqrt(2), atol=1e-10)
    @test z == 0.0
    
    # Test magnetic field magnitude
    B_mag = magnetic_field_magnitude(1.0, 2.0, 3.0)
    @test isapprox(B_mag, sqrt(14), atol=1e-10)
    
    # Test safety factor
    q = safety_factor(1.0, 2.0, 1.0)
    @test q == 2.0
    
    println("✓ Utility functions test passed")
end

println("\nAll tests passed! ✓")
println("StellaratorPhysics package is working correctly.")
