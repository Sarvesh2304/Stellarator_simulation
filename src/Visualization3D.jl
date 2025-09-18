"""
3D Visualization Module for Stellarator Physics

This module provides comprehensive 3D visualization capabilities for:
- Magnetic field lines and surfaces
- Plasma density and temperature profiles
- Transport coefficient visualization
- Optimization results visualization
"""

using Plots
using PlotlyJS
using LinearAlgebra

export Visualization3D, plot_magnetic_field, plot_plasma_surfaces
export plot_transport_coefficients, plot_optimization_results
export plot_field_lines, plot_magnetic_surfaces, plot_3d_geometry

"""
    struct Visualization3D

Main structure for 3D visualization capabilities.
"""
mutable struct Visualization3D
    bfield::MagneticField3D
    transport::NeoclassicalTransport
    plot_style::String
    color_scheme::String
    
    function Visualization3D(bfield::MagneticField3D, transport::NeoclassicalTransport)
        new(bfield, transport, "plotly", "viridis")
    end
end

"""
    plot_magnetic_field(bfield::MagneticField3D, R_range, φ_range, Z_range; n_points=50)

Plot 3D magnetic field vectors and magnitude.
"""
function plot_magnetic_field(bfield::MagneticField3D, R_range, φ_range, Z_range; n_points=50)
    # Create grid
    R_vals = range(R_range[1], R_range[2], length=n_points)
    φ_vals = range(φ_range[1], φ_range[2], length=n_points)
    Z_vals = range(Z_range[1], Z_range[2], length=n_points)
    
    # Calculate magnetic field on grid
    B_magnitude = zeros(n_points, n_points, n_points)
    B_x = zeros(n_points, n_points, n_points)
    B_y = zeros(n_points, n_points, n_points)
    B_z = zeros(n_points, n_points, n_points)
    
    for i in 1:n_points
        for j in 1:n_points
            for k in 1:n_points
                R = R_vals[i]
                φ = φ_vals[j]
                Z = Z_vals[k]
                
                B_R, B_φ, B_Z = calculate_magnetic_field(bfield, R, φ, Z)
                B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
                
                B_magnitude[i, j, k] = B_mag
                
                # Convert to Cartesian
                B_x[i, j, k] = B_R * cos(φ) - B_φ * sin(φ)
                B_y[i, j, k] = B_R * sin(φ) + B_φ * cos(φ)
                B_z[i, j, k] = B_Z
            end
        end
    end
    
    # Create 3D plot
    plot = plot3d(
        xlabel="X [m]",
        ylabel="Y [m]",
        zlabel="Z [m]",
        title="3D Magnetic Field Magnitude",
        colorbar_title="|B| [T]"
    )
    
    # Add isosurface
    add_trace!(plot, isosurface(
        x=R_vals,
        y=φ_vals,
        z=Z_vals,
        value=B_magnitude,
        isomin=minimum(B_magnitude),
        isomax=maximum(B_magnitude),
        surface_count=10,
        colorscale="viridis"
    ))
    
    return plot
end

"""
    plot_plasma_surfaces(bfield::MagneticField3D, s_values; n_θ=64, n_φ=64)

Plot multiple magnetic surfaces.
"""
function plot_plasma_surfaces(bfield::MagneticField3D, s_values; n_θ=64, n_φ=64)
    plot = plot3d(
        xlabel="X [m]",
        ylabel="Y [m]",
        zlabel="Z [m]",
        title="Magnetic Surfaces",
        showlegend=true
    )
    
    colors = ["red", "blue", "green", "orange", "purple", "brown", "pink", "gray"]
    
    for (i, s) in enumerate(s_values)
        surface = find_magnetic_surface(bfield, s, n_θ, n_φ)
        
        # Extract coordinates
        x_vals = [point[1] for point in surface.points]
        y_vals = [point[2] for point in surface.points]
        z_vals = [point[3] for point in surface.points]
        
        # Create surface plot
        add_trace!(plot, mesh3d(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            opacity=0.3,
            color=colors[mod(i-1, length(colors)) + 1],
            name="s = $(round(s, digits=2))"
        ))
    end
    
    return plot
end

"""
    plot_field_lines(bfield::MagneticField3D, start_points; length=10.0, n_steps=1000)

Plot magnetic field lines starting from given points.
"""
function plot_field_lines(bfield::MagneticField3D, start_points; length=10.0, n_steps=1000)
    plot = plot3d(
        xlabel="X [m]",
        ylabel="Y [m]",
        zlabel="Z [m]",
        title="Magnetic Field Lines",
        showlegend=true
    )
    
    colors = ["red", "blue", "green", "orange", "purple", "brown", "pink", "gray"]
    
    for (i, start_point) in enumerate(start_points)
        R₀, φ₀, Z₀ = start_point
        
        # Trace field line
        field_line = trace_field_line(bfield, R₀, φ₀, Z₀, length)
        
        # Extract coordinates
        x_vals = [point[1] for point in field_line]
        y_vals = [point[2] for point in field_line]
        z_vals = [point[3] for point in field_line]
        
        # Add field line
        add_trace!(plot, scatter3d(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            mode="lines",
            line=attr(color=colors[mod(i-1, length(colors)) + 1], width=3),
            name="Field Line $(i)"
        ))
    end
    
    return plot
end

"""
    plot_transport_coefficients(transport::NeoclassicalTransport, s_values)

Plot transport coefficients as a function of radius.
"""
function plot_transport_coefficients(transport::NeoclassicalTransport, s_values)
    # Calculate transport coefficients
    D_11_vals = Float64[]
    D_22_vals = Float64[]
    χ_eff_vals = Float64[]
    j_bs_vals = Float64[]
    
    for s in s_values
        coeffs = calculate_transport_coefficients(transport, s)
        push!(D_11_vals, coeffs.D_11)
        push!(D_22_vals, coeffs.D_22)
        push!(χ_eff_vals, coeffs.χ_eff)
        
        j_bs = bootstrap_current(transport, s)
        push!(j_bs_vals, j_bs)
    end
    
    # Create subplots
    plot1 = plot(
        s_values, D_11_vals,
        xlabel="Normalized Radius (s)",
        ylabel="D₁₁ [m²/s]",
        title="Particle Diffusion Coefficient",
        linewidth=2,
        color=:blue
    )
    
    plot2 = plot(
        s_values, D_22_vals,
        xlabel="Normalized Radius (s)",
        ylabel="D₂₂ [m²/s]",
        title="Heat Diffusion Coefficient",
        linewidth=2,
        color=:red
    )
    
    plot3 = plot(
        s_values, χ_eff_vals,
        xlabel="Normalized Radius (s)",
        ylabel="χ_eff [m²/s]",
        title="Effective Thermal Diffusivity",
        linewidth=2,
        color=:green
    )
    
    plot4 = plot(
        s_values, j_bs_vals,
        xlabel="Normalized Radius (s)",
        ylabel="j_bs [A/m²]",
        title="Bootstrap Current Density",
        linewidth=2,
        color=:orange
    )
    
    # Combine plots
    combined_plot = plot(plot1, plot2, plot3, plot4, layout=(2, 2), size=(800, 600))
    
    return combined_plot
end

"""
    plot_optimization_results(bfield::MagneticField3D, optimization_results)

Plot optimization results and convergence.
"""
function plot_optimization_results(bfield::MagneticField3D, optimization_results)
    # Plot objective function convergence
    plot1 = plot(
        optimization_results.objective_values,
        xlabel="Iteration",
        ylabel="Objective Function Value",
        title="Optimization Convergence",
        linewidth=2,
        color=:blue
    )
    
    # Plot parameter evolution
    plot2 = plot(
        optimization_results.parameter_history,
        xlabel="Iteration",
        ylabel="Parameter Value",
        title="Parameter Evolution",
        linewidth=2,
        color=:red
    )
    
    # Plot magnetic field before and after optimization
    plot3 = plot_magnetic_field(bfield, (0.5, 1.5), (0, 2π), (-0.5, 0.5))
    
    # Combine plots
    combined_plot = plot(plot1, plot2, plot3, layout=(2, 2), size=(1000, 800))
    
    return combined_plot
end

"""
    plot_3d_geometry(bfield::MagneticField3D, s_values; n_θ=64, n_φ=64)

Plot 3D stellarator geometry with magnetic surfaces.
"""
function plot_3d_geometry(bfield::MagneticField3D, s_values; n_θ=64, n_φ=64)
    plot = plot3d(
        xlabel="X [m]",
        ylabel="Y [m]",
        zlabel="Z [m]",
        title="3D Stellarator Geometry",
        showlegend=true,
        aspect_ratio=:equal
    )
    
    # Plot magnetic surfaces
    for (i, s) in enumerate(s_values)
        surface = find_magnetic_surface(bfield, s, n_θ, n_φ)
        
        # Extract coordinates
        x_vals = [point[1] for point in surface.points]
        y_vals = [point[2] for point in surface.points]
        z_vals = [point[3] for point in surface.points]
        
        # Create surface plot
        add_trace!(plot, mesh3d(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            opacity=0.3,
            color="rgba(0, 0, 255, 0.3)",
            name="Magnetic Surface s=$(round(s, digits=2))"
        ))
    end
    
    # Add coordinate axes
    add_trace!(plot, scatter3d(
        x=[0, 2], y=[0, 0], z=[0, 0],
        mode="lines",
        line=attr(color="red", width=5),
        name="X-axis"
    ))
    
    add_trace!(plot, scatter3d(
        x=[0, 0], y=[0, 2], z=[0, 0],
        mode="lines",
        line=attr(color="green", width=5),
        name="Y-axis"
    ))
    
    add_trace!(plot, scatter3d(
        x=[0, 0], y=[0, 0], z=[0, 2],
        mode="lines",
        line=attr(color="blue", width=5),
        name="Z-axis"
    ))
    
    return plot
end

"""
    plot_quasi_symmetry_measure(bfield::MagneticField3D, s_values)

Plot quasi-symmetry measure as a function of radius.
"""
function plot_quasi_symmetry_measure(bfield::MagneticField3D, s_values)
    qs_measures = Float64[]
    
    for s in s_values
        surface = find_magnetic_surface(bfield, s)
        measure = quasi_symmetry_measure(bfield, surface)
        push!(qs_measures, measure)
    end
    
    plot = plot(
        s_values, qs_measures,
        xlabel="Normalized Radius (s)",
        ylabel="Quasi-Symmetry Measure",
        title="Quasi-Symmetry Measure vs Radius",
        linewidth=2,
        color=:purple
    )
    
    return plot
end

"""
    plot_magnetic_well(bfield::MagneticField3D, s_values)

Plot magnetic well depth as a function of radius.
"""
function plot_magnetic_well(bfield::MagneticField3D, s_values)
    well_depths = Float64[]
    
    for s in s_values
        well_depth = magnetic_well_depth(bfield, 0.1, s)
        push!(well_depths, well_depth)
    end
    
    plot = plot(
        s_values, well_depths,
        xlabel="Normalized Radius (s)",
        ylabel="Magnetic Well Depth",
        title="Magnetic Well Depth vs Radius",
        linewidth=2,
        color=:brown
    )
    
    return plot
end

"""
    plot_comparison_results(comparison_results)

Plot comparison results between stellarator and tokamak.
"""
function plot_comparison_results(comparison_results)
    # Extract data for plotting
    metrics = String[]
    stellarator_values = Float64[]
    tokamak_values = Float64[]
    ratios = Float64[]
    
    for (category, results) in comparison_results
        for result in results
            push!(metrics, result.metric)
            push!(stellarator_values, result.stellarator_value)
            push!(tokamak_values, result.tokamak_value)
            push!(ratios, result.ratio)
        end
    end
    
    # Create comparison plot
    plot1 = bar(
        metrics, ratios,
        xlabel="Metric",
        ylabel="Stellarator/Tokamak Ratio",
        title="Performance Comparison",
        color=:blue,
        xtickangle=45
    )
    
    # Add horizontal line at ratio = 1
    hline!(plot1, [1.0], color=:red, linestyle=:dash, linewidth=2)
    
    # Create scatter plot
    plot2 = scatter(
        tokamak_values, stellarator_values,
        xlabel="Tokamak Value",
        ylabel="Stellarator Value",
        title="Stellarator vs Tokamak Values",
        color=:green,
        markersize=6
    )
    
    # Add diagonal line
    min_val = min(minimum(tokamak_values), minimum(stellarator_values))
    max_val = max(maximum(tokamak_values), maximum(stellarator_values))
    plot!(plot2, [min_val, max_val], [min_val, max_val], 
          color=:red, linestyle=:dash, linewidth=2)
    
    # Combine plots
    combined_plot = plot(plot1, plot2, layout=(1, 2), size=(1200, 600))
    
    return combined_plot
end
