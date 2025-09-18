"""
Comprehensive Stellarator Physics Analysis Example

This example demonstrates the complete functionality of the StellaratorPhysics package,
including:
- 3D magnetic field calculations
- Neoclassical transport analysis
- Stellarator optimization
- Comparison with tokamak performance
- 3D visualization
"""

using Pkg
Pkg.activate(".")

using StellaratorPhysics
using Plots
using PlotlyJS

# Set up plotting backend
plotlyjs()

println("=== Stellarator Physics Analysis Example ===")
println()

# 1. Create a stellarator magnetic field configuration
println("1. Creating stellarator magnetic field configuration...")
R₀ = 1.0  # Major radius [m]
a = 0.2   # Minor radius [m]
N = 5     # Number of field periods
B₀ = 1.0  # Reference magnetic field [T]

stellarator_bfield = MagneticField3D(R₀, a, N, B₀)

# Add stellarator harmonics
stellarator_bfield.harmonics = stellarator_harmonics(5, 5, R₀, a, N)

println("   ✓ Stellarator field created with R₀ = $(R₀) m, a = $(a) m, N = $(N)")

# 2. Create a tokamak magnetic field for comparison
println("\n2. Creating tokamak magnetic field for comparison...")
tokamak_bfield = create_tokamak_field(R₀, a, 1.0, 2.0, B₀)
println("   ✓ Tokamak field created for comparison")

# 3. Set up plasma parameters
println("\n3. Setting up plasma parameters...")
T_e = 1000.0  # Electron temperature [eV]
T_i = 1000.0  # Ion temperature [eV]
n_e = 1e20    # Electron density [m^-3]
n_i = 1e20    # Ion density [m^-3]
Z_eff = 1.0   # Effective charge

stellarator_transport = NeoclassicalTransport(stellarator_bfield, T_e, T_i, n_e, n_i, Z_eff)
tokamak_transport = NeoclassicalTransport(tokamak_bfield, T_e, T_i, n_e, n_i, Z_eff)

println("   ✓ Plasma parameters set: T_e = $(T_e) eV, n_e = $(n_e) m^-3")

# 4. Calculate magnetic field properties
println("\n4. Calculating magnetic field properties...")
s_values = range(0.1, 0.9, length=20)

# Calculate safety factors
q_stellarator = Float64[]
q_tokamak = Float64[]

for s in s_values
    surface_s = find_magnetic_surface(stellarator_bfield, s)
    surface_t = find_magnetic_surface(tokamak_bfield, s)
    
    push!(q_stellarator, surface_s.q)
    push!(q_tokamak, surface_t.q)
end

println("   ✓ Safety factors calculated for both configurations")

# 5. Calculate transport coefficients
println("\n5. Calculating neoclassical transport coefficients...")
D_11_stellarator = Float64[]
D_22_stellarator = Float64[]
χ_eff_stellarator = Float64[]

D_11_tokamak = Float64[]
D_22_tokamak = Float64[]
χ_eff_tokamak = Float64[]

for s in s_values
    coeffs_s = calculate_transport_coefficients(stellarator_transport, s)
    coeffs_t = calculate_transport_coefficients(tokamak_transport, s)
    
    push!(D_11_stellarator, coeffs_s.D_11)
    push!(D_22_stellarator, coeffs_s.D_22)
    push!(χ_eff_stellarator, coeffs_s.χ_eff)
    
    push!(D_11_tokamak, coeffs_t.D_11)
    push!(D_22_tokamak, coeffs_t.D_22)
    push!(χ_eff_tokamak, coeffs_t.χ_eff)
end

println("   ✓ Transport coefficients calculated")

# 6. Perform stellarator optimization
println("\n6. Performing stellarator optimization...")
println("   Optimizing for quasi-symmetry...")

try
    result, optimal_harmonics = optimize_quasi_symmetry(stellarator_bfield, 100)
    
    if result.success
        println("   ✓ Optimization completed successfully")
        println("   Final objective value: $(result.objective_value)")
    else
        println("   ⚠ Optimization did not converge")
    end
catch e
    println("   ⚠ Optimization failed: $(e)")
end

# 7. Compare stellarator and tokamak performance
println("\n7. Comparing stellarator and tokamak performance...")
comparison = TokamakComparison(stellarator_bfield, tokamak_bfield, 
                              stellarator_transport, tokamak_transport)

# Compare at mid-radius
s_compare = 0.5
performance_results = compare_performance(comparison, s_compare)
transport_results = compare_transport(comparison, s_compare)
stability_results = compare_stability(comparison, s_compare)

println("   ✓ Performance comparison completed")

# 8. Create visualizations
println("\n8. Creating visualizations...")

# Plot safety factor profiles
plot1 = plot(
    s_values, q_stellarator,
    label="Stellarator",
    xlabel="Normalized Radius (s)",
    ylabel="Safety Factor (q)",
    title="Safety Factor Profiles",
    linewidth=2,
    color=:blue
)
plot!(plot1, s_values, q_tokamak, label="Tokamak", linewidth=2, color=:red)

# Plot transport coefficients
plot2 = plot(
    s_values, D_11_stellarator,
    label="Stellarator D₁₁",
    xlabel="Normalized Radius (s)",
    ylabel="D₁₁ [m²/s]",
    title="Particle Diffusion Coefficient",
    linewidth=2,
    color=:blue
)
plot!(plot2, s_values, D_11_tokamak, label="Tokamak D₁₁", linewidth=2, color=:red)

plot3 = plot(
    s_values, χ_eff_stellarator,
    label="Stellarator χ_eff",
    xlabel="Normalized Radius (s)",
    ylabel="χ_eff [m²/s]",
    title="Effective Thermal Diffusivity",
    linewidth=2,
    color=:blue
)
plot!(plot3, s_values, χ_eff_tokamak, label="Tokamak χ_eff", linewidth=2, color=:red)

# Combine plots
combined_plot = plot(plot1, plot2, plot3, layout=(1, 3), size=(1200, 400))

# Save plots
savefig(combined_plot, "stellarator_analysis_results.png")
println("   ✓ Plots saved as 'stellarator_analysis_results.png'")

# 9. Create 3D visualizations
println("\n9. Creating 3D visualizations...")

# Plot magnetic surfaces
s_surfaces = [0.2, 0.4, 0.6, 0.8]
plot_3d = plot_plasma_surfaces(stellarator_bfield, s_surfaces)

# Save 3D plot
savefig(plot_3d, "stellarator_3d_surfaces.html")
println("   ✓ 3D magnetic surfaces saved as 'stellarator_3d_surfaces.html'")

# 10. Print summary results
println("\n=== ANALYSIS SUMMARY ===")
println()

# Safety factor comparison
q_s_mid = q_stellarator[Int(length(q_stellarator)/2)]
q_t_mid = q_tokamak[Int(length(q_tokamak)/2)]
println("Safety Factor at mid-radius:")
println("  Stellarator: $(round(q_s_mid, digits=3))")
println("  Tokamak: $(round(q_t_mid, digits=3))")
println("  Ratio: $(round(q_s_mid/q_t_mid, digits=3))")

# Transport coefficient comparison
D_11_s_mid = D_11_stellarator[Int(length(D_11_stellarator)/2)]
D_11_t_mid = D_11_tokamak[Int(length(D_11_tokamak)/2)]
println("\nParticle Diffusion Coefficient at mid-radius:")
println("  Stellarator: $(round(D_11_s_mid, digits=6)) m²/s")
println("  Tokamak: $(round(D_11_t_mid, digits=6)) m²/s")
println("  Ratio: $(round(D_11_s_mid/D_11_t_mid, digits=3))")

χ_eff_s_mid = χ_eff_stellarator[Int(length(χ_eff_stellarator)/2)]
χ_eff_t_mid = χ_eff_tokamak[Int(length(χ_eff_tokamak)/2)]
println("\nThermal Diffusivity at mid-radius:")
println("  Stellarator: $(round(χ_eff_s_mid, digits=6)) m²/s")
println("  Tokamak: $(round(χ_eff_t_mid, digits=6)) m²/s")
println("  Ratio: $(round(χ_eff_s_mid/χ_eff_t_mid, digits=3))")

# Performance comparison results
println("\nPerformance Comparison at mid-radius:")
for result in performance_results
    println("  $(result.metric):")
    println("    Stellarator: $(round(result.stellarator_value, digits=6))")
    println("    Tokamak: $(round(result.tokamak_value, digits=6))")
    println("    Ratio: $(round(result.ratio, digits=3)) ($(result.advantage) advantage)")
end

println("\n=== ANALYSIS COMPLETE ===")
println("Results saved:")
println("  - stellarator_analysis_results.png (2D plots)")
println("  - stellarator_3d_surfaces.html (3D visualization)")
println()
println("This example demonstrates the complete functionality of the StellaratorPhysics package.")
println("You can modify the parameters above to explore different configurations.")
