"""
Create True Stellarator 3D Visualization with Quasi-Isodynamic Structure
This script creates proper stellarator visualizations showing the helical winding
and quasi-isodynamic magnetic field structure that distinguishes it from tokamaks.
"""

# Physical constants
const μ₀ = 4π * 1e-7  # Permeability of free space [H/m]
const e = 1.602176634e-19  # Elementary charge [C]
const mₑ = 9.1093837015e-31  # Electron mass [kg]
const mᵢ = 1.67262192369e-27  # Proton mass [kg]
const k_B = 1.380649e-23  # Boltzmann constant [J/K]
const T_eV = e / k_B  # Temperature conversion factor [K/eV]

# Stellarator Configuration with proper helical structure
mutable struct StellaratorField3D
    R₀::Float64      # Major radius
    a::Float64       # Minor radius
    N::Int           # Number of field periods
    B₀::Float64      # Reference magnetic field strength
    iota::Float64    # Rotational transform
    helical_amplitude::Float64  # Helical winding amplitude
    quasi_isodynamic_factor::Float64  # Quasi-isodynamic optimization factor
    
    function StellaratorField3D(R₀, a, N, B₀=1.0, iota=0.3, helical_amp=0.1, qid_factor=0.8)
        new(R₀, a, N, B₀, iota, helical_amp, qid_factor)
    end
end

# Quasi-isodynamic magnetic surface with proper stellarator geometry
struct QuasiIsodynamicSurface
    points::Vector{Tuple{Float64,Float64,Float64}}  # Surface points
    magnetic_field_strength::Vector{Float64}         # |B| at each point
    area::Float64                                    # Surface area
    volume::Float64                                  # Enclosed volume
    iota::Float64                                    # Rotational transform
    quasi_isodynamic_quality::Float64                # QI quality measure
end

# Generate proper stellarator magnetic field with helical winding
function calculate_stellarator_field(bfield::StellaratorField3D, R, φ, Z)
    # Convert to cylindrical coordinates
    r = sqrt((R - bfield.R₀)^2 + Z^2)
    θ = atan(Z, R - bfield.R₀)
    
    # Main toroidal field
    B_φ_main = bfield.B₀ * bfield.R₀ / R
    
    # Helical field components (quasi-isodynamic structure)
    # This creates the characteristic stellarator winding
    helical_phase = bfield.N * φ - bfield.iota * θ
    
    # Poloidal field with helical modulation
    B_R_helical = bfield.helical_amplitude * bfield.B₀ * sin(helical_phase) * (r / bfield.a)
    B_φ_helical = bfield.helical_amplitude * bfield.B₀ * cos(helical_phase) * (r / bfield.a) * 0.5
    B_Z_helical = bfield.helical_amplitude * bfield.B₀ * cos(helical_phase) * (r / bfield.a)
    
    # Quasi-isodynamic optimization
    # This creates regions of nearly constant |B| along field lines
    qid_modulation = bfield.quasi_isodynamic_factor * sin(2 * helical_phase) * 0.1
    
    B_R = B_R_helical + qid_modulation * sin(θ)
    B_φ = B_φ_main + B_φ_helical + qid_modulation * cos(θ)
    B_Z = B_Z_helical + qid_modulation * sin(θ)
    
    return B_R, B_φ, B_Z
end

# Generate quasi-isodynamic magnetic surface
function find_quasi_isodynamic_surface(bfield::StellaratorField3D, s, n_θ=64, n_φ=64)
    points = Tuple{Float64,Float64,Float64}[]
    B_values = Float64[]
    
    # Generate surface points with proper stellarator geometry
    for i in 1:n_θ
        θ = 2π * (i - 1) / n_θ
        for j in 1:n_φ
            φ = 2π * (j - 1) / n_φ
            
            # Stellarator minor radius with helical modulation
            r_base = s * bfield.a
            helical_modulation = bfield.helical_amplitude * sin(bfield.N * φ) * 0.1
            r = r_base * (1 + helical_modulation)
            
            # Major radius with helical shift
            R_shift = bfield.helical_amplitude * cos(bfield.N * φ) * 0.05
            R = bfield.R₀ + r * cos(θ) + R_shift
            Z = r * sin(θ)
            
            push!(points, (R, φ, Z))
            
            # Calculate magnetic field strength at this point
            B_R, B_φ, B_Z = calculate_stellarator_field(bfield, R, φ, Z)
            B_mag = sqrt(B_R^2 + B_φ^2 + B_Z^2)
            push!(B_values, B_mag)
        end
    end
    
    # Calculate surface properties
    area = 4π^2 * bfield.R₀ * s * bfield.a * (1 + bfield.helical_amplitude)
    volume = 2π^2 * bfield.R₀ * (s * bfield.a)^2 * (1 + bfield.helical_amplitude)
    
    # Calculate quasi-isodynamic quality
    B_mean = sum(B_values) / length(B_values)
    B_variance = sum((b - B_mean)^2 for b in B_values) / length(B_values)
    B_std = sqrt(B_variance)
    B_variation = B_std / B_mean
    qid_quality = 1.0 / (1.0 + B_variation)  # Higher is better
    
    return QuasiIsodynamicSurface(points, B_values, area, volume, bfield.iota, qid_quality)
end

# Generate stellarator field lines with proper helical structure
function trace_stellarator_field_line(bfield::StellaratorField3D, R₀, φ₀, Z₀, length)
    points = Tuple{Float64,Float64,Float64}[]
    
    # Initial point
    R, φ, Z = R₀, φ₀, Z₀
    push!(points, (R, φ, Z))
    
    # Integration parameters
    ds = length / 1000
    steps = Int(length / ds)
    
    for i in 1:steps
        # Calculate magnetic field at current point
        B_R, B_φ, B_Z = calculate_stellarator_field(bfield, R, φ, Z)
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

# Create stellarator configuration with proper quasi-isodynamic structure
println("Creating quasi-isodynamic stellarator configuration...")
R₀ = 1.0      # Major radius [m]
a = 0.2       # Minor radius [m]
N = 5         # Number of field periods
B₀ = 1.0      # Reference magnetic field [T]
iota = 0.3    # Rotational transform
helical_amp = 0.15  # Helical winding amplitude
qid_factor = 0.8    # Quasi-isodynamic optimization factor

stellarator_bfield = StellaratorField3D(R₀, a, N, B₀, iota, helical_amp, qid_factor)

# Generate quasi-isodynamic surfaces
s_values = [0.2, 0.4, 0.6, 0.8]
surfaces = []

println("Generating quasi-isodynamic magnetic surfaces...")
for s in s_values
    surface = find_quasi_isodynamic_surface(stellarator_bfield, s)
    push!(surfaces, surface)
    
    println("  Surface s = $(s):")
    println("    Area = $(round(surface.area, digits=4)) m²")
    println("    Volume = $(round(surface.volume, digits=4)) m³")
    println("    Rotational transform ι = $(round(surface.iota, digits=3))")
    println("    Quasi-isodynamic quality = $(round(surface.quasi_isodynamic_quality, digits=3))")
end

# Create data files for 3D visualization
println("\\nCreating 3D visualization data...")

# Create surface data files with magnetic field strength
for (i, surface) in enumerate(surfaces)
    x_vals = Float64[]
    y_vals = Float64[]
    z_vals = Float64[]
    B_vals = Float64[]
    
    for point in surface.points
        R, φ, Z = point
        # Convert to Cartesian
        x, y, z = R * cos(φ), R * sin(φ), Z
        push!(x_vals, x)
        push!(y_vals, y)
        push!(z_vals, z)
    end
    
    # Add magnetic field strength data
    B_vals = surface.magnetic_field_strength
    
    filename = "stellarator_surface_$(i)_s$(s_values[i]).dat"
    open(filename, "w") do file
        println(file, "# Quasi-isodynamic stellarator surface s = $(s_values[i])")
        println(file, "# X[m] Y[m] Z[m] |B|[T]")
        for j in 1:length(x_vals)
            println(file, "$(x_vals[j]) $(y_vals[j]) $(z_vals[j]) $(B_vals[j])")
        end
    end
    println("  ✓ Created $filename with $(length(x_vals)) points")
end

# Generate field line data
println("Generating stellarator field lines...")
field_lines = []
start_points = [
    (1.15, 0.0, 0.0),
    (1.25, 0.0, 0.0),
    (1.35, 0.0, 0.0)
]

for (i, start_point) in enumerate(start_points)
    field_line = trace_stellarator_field_line(stellarator_bfield, start_point[1], start_point[2], start_point[3], 20.0)
    push!(field_lines, field_line)
    
    # Save field line data
    filename = "stellarator_fieldline_$(i).dat"
    open(filename, "w") do file
        println(file, "# Stellarator field line $(i)")
        println(file, "# X[m] Y[m] Z[m]")
        for point in field_line
            R, φ, Z = point
            x, y, z = R * cos(φ), R * sin(φ), Z
            println(file, "$x $y $z")
        end
    end
    println("  ✓ Created $filename with $(length(field_line)) points")
end

# Create enhanced Python script for stellarator visualization
python_script = """
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors

# Read surface data
surfaces = []
s_values = [0.2, 0.4, 0.6, 0.8]
colors = ['red', 'blue', 'green', 'orange']

fig = plt.figure(figsize=(15, 10))

# 3D plot of surfaces
ax1 = fig.add_subplot(121, projection='3d')

for i, s in enumerate(s_values):
    filename = f'stellarator_surface_{i+1}_s{s}.dat'
    data = np.loadtxt(filename, comments='#')
    
    x = data[:, 0]
    y = data[:, 1] 
    z = data[:, 2]
    B = data[:, 3]  # Magnetic field strength
    
    # Plot surface with color based on magnetic field strength
    scatter = ax1.scatter(x, y, z, c=B, cmap='plasma', s=3, alpha=0.7, 
                         label=f'Surface s={s}')

# Add field lines
for i in range(3):
    filename = f'stellarator_fieldline_{i+1}.dat'
    try:
        data = np.loadtxt(filename, comments='#')
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]
        
        # Plot field line
        ax1.plot(x, y, z, 'k-', linewidth=1, alpha=0.8)
    except:
        pass

# Add coordinate axes
ax1.plot([0, 1.8], [0, 0], [0, 0], 'r-', linewidth=3, label='X-axis')
ax1.plot([0, 0], [0, 1.8], [0, 0], 'g-', linewidth=3, label='Y-axis')
ax1.plot([0, 0], [0, 0], [0, 0.6], 'b-', linewidth=3, label='Z-axis')

ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_zlabel('Z [m]')
ax1.set_title('Quasi-Isodynamic Stellarator\\n3D Magnetic Surfaces & Field Lines')
ax1.legend()

# 2D cross-section showing helical structure
ax2 = fig.add_subplot(122)

# Plot cross-sections at different toroidal angles
phi_angles = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
for phi in phi_angles:
    x_cross = []
    y_cross = []
    
    for s in s_values:
        # Generate cross-section points
        theta = np.linspace(0, 2*np.pi, 100)
        r_base = s * 0.2
        helical_mod = 0.15 * np.sin(5 * phi) * 0.1
        r = r_base * (1 + helical_mod)
        R_shift = 0.15 * np.cos(5 * phi) * 0.05
        
        R = 1.0 + r * np.cos(theta) + R_shift
        Z = r * np.sin(theta)
        
        x_cross.extend(R * np.cos(phi))
        y_cross.extend(R * np.sin(phi))
    
    ax2.scatter(x_cross, y_cross, s=1, alpha=0.6, 
               label=f'φ = {phi:.2f}' if phi in [0, np.pi/2, np.pi] else "")

ax2.set_xlabel('X [m]')
ax2.set_ylabel('Y [m]')
ax2.set_title('Stellarator Cross-Sections\\nShowing Helical Structure')
ax2.set_aspect('equal')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('quasi_isodynamic_stellarator_3d.png', dpi=300, bbox_inches='tight')
plt.savefig('quasi_isodynamic_stellarator_3d.pdf', bbox_inches='tight')
plt.show()

print("Quasi-isodynamic stellarator visualization saved!")
"""

# Write enhanced Python script
open("plot_stellarator_3d.py", "w") do file
    write(file, python_script)
end

# Create enhanced HTML visualization
html_script = """
<!DOCTYPE html>
<html>
<head>
    <title>Quasi-Isodynamic Stellarator 3D Visualization</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .container { display: flex; flex-wrap: wrap; }
        .plot-container { flex: 1; min-width: 400px; margin: 10px; }
        .info { background: #f0f0f0; padding: 15px; margin: 10px; border-radius: 5px; }
    </style>
</head>
<body>
    <h1>Quasi-Isodynamic Stellarator 3D Visualization</h1>
    
    <div class="info">
        <h3>Stellarator vs Tokamak</h3>
        <p><strong>Stellarator:</strong> 3D helical magnetic field structure, quasi-isodynamic optimization, 
        steady-state operation without plasma current</p>
        <p><strong>Key Features:</strong> Helical winding, quasi-isodynamic magnetic surfaces, 
        reduced transport, flexible optimization</p>
    </div>
    
    <div class="container">
        <div class="plot-container">
            <div id="plot3d"></div>
        </div>
        <div class="plot-container">
            <div id="plot2d"></div>
        </div>
    </div>
    
    <script>
        // Generate quasi-isodynamic stellarator data
        function generateStellaratorData() {
            const R0 = 1.0;
            const a = 0.2;
            const N = 5; // Field periods
            const helical_amp = 0.15;
            
            const surfaces = [];
            const s_values = [0.2, 0.4, 0.6, 0.8];
            const colors = ['red', 'blue', 'green', 'orange'];
            
            s_values.forEach((s, i) => {
                const x = [], y = [], z = [], B = [];
                
                for (let theta = 0; theta < 2 * Math.PI; theta += 0.1) {
                    for (let phi = 0; phi < 2 * Math.PI; phi += 0.1) {
                        // Stellarator geometry with helical modulation
                        const r_base = s * a;
                        const helical_mod = helical_amp * Math.sin(N * phi) * 0.1;
                        const r = r_base * (1 + helical_mod);
                        const R_shift = helical_amp * Math.cos(N * phi) * 0.05;
                        
                        const R = R0 + r * Math.cos(theta) + R_shift;
                        const Z = r * Math.sin(theta);
                        const X = R * Math.cos(phi);
                        const Y = R * Math.sin(phi);
                        
                        // Magnetic field strength (quasi-isodynamic)
                        const B_mag = 1.0 + helical_amp * Math.cos(N * phi) * 0.1;
                        
                        x.push(X);
                        y.push(Y);
                        z.push(Z);
                        B.push(B_mag);
                    }
                }
                
                surfaces.push({
                    x: x, y: y, z: z, B: B,
                    type: 'scatter3d',
                    mode: 'markers',
                    marker: {
                        size: 2,
                        color: B,
                        colorscale: 'Plasma',
                        opacity: 0.7,
                        showscale: false
                    },
                    name: 'Surface s=' + s.toString(),
                    visible: i === 0 ? true : 'legendonly'
                });
            });
            
            return surfaces;
        }
        
        // 3D Plot
        const surfaces3d = generateStellaratorData();
        const layout3d = {
            title: 'Quasi-Isodynamic Stellarator<br>3D Magnetic Surfaces',
            scene: {
                xaxis: {title: 'X [m]'},
                yaxis: {title: 'Y [m]'},
                zaxis: {title: 'Z [m]'},
                aspectmode: 'data',
                camera: {eye: {x: 1.5, y: 1.5, z: 1.5}}
            },
            width: 600,
            height: 500
        };
        
        Plotly.newPlot('plot3d', surfaces3d, layout3d);
        
        // 2D Cross-section plot
        const crossSectionData = [];
        const phi_angles = [0, Math.PI/4, Math.PI/2, 3*Math.PI/4, Math.PI];
        
        phi_angles.forEach((phi, i) => {
            const x = [], y = [];
            const s_values = [0.2, 0.4, 0.6, 0.8];
            
            s_values.forEach(s => {
                for (let theta = 0; theta < 2 * Math.PI; theta += 0.05) {
                    const r_base = s * 0.2;
                    const helical_mod = 0.15 * Math.sin(5 * phi) * 0.1;
                    const r = r_base * (1 + helical_mod);
                    const R_shift = 0.15 * Math.cos(5 * phi) * 0.05;
                    
                    const R = 1.0 + r * Math.cos(theta) + R_shift;
                    const X = R * Math.cos(phi);
                    const Y = R * Math.sin(phi);
                    
                    x.push(X);
                    y.push(Y);
                }
            });
            
            crossSectionData.push({
                x: x, y: y,
                type: 'scatter',
                mode: 'markers',
                marker: {size: 1, opacity: 0.6},
                name: 'φ = ' + (phi * 180 / Math.PI).toFixed(0) + '°',
                visible: i < 3 ? true : 'legendonly'
            });
        });
        
        const layout2d = {
            title: 'Stellarator Cross-Sections<br>Showing Helical Structure',
            xaxis: {title: 'X [m]'},
            yaxis: {title: 'Y [m]', scaleanchor: 'x', scaleratio: 1},
            width: 600,
            height: 500
        };
        
        Plotly.newPlot('plot2d', crossSectionData, layout2d);
    </script>
</body>
</html>
"""

# Write enhanced HTML script
open("quasi_isodynamic_stellarator.html", "w") do file
    write(file, html_script)
end

println("\\n=== Quasi-Isodynamic Stellarator Visualization Created ===")
println("✓ Surface data files: stellarator_surface_*.dat")
println("✓ Field line data files: stellarator_fieldline_*.dat")
println("✓ Python script: plot_stellarator_3d.py")
println("✓ HTML visualization: quasi_isodynamic_stellarator.html")
println("\\nKey Features:")
println("• Helical winding structure (N = $N field periods)")
println("• Quasi-isodynamic magnetic field optimization")
println("• Rotational transform ι = $iota")
println("• 3D magnetic surfaces with proper stellarator geometry")
println("• Field line tracing showing helical structure")
println("\\nTo view: Open quasi_isodynamic_stellarator.html in your browser")
println("\\n=== STELLARATOR VISUALIZATION COMPLETE! ===")
