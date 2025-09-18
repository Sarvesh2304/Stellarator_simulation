import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors

# Read surface data
surfaces = []
s_values = [0.2, 0.4, 0.6, 0.8]
colors = ['red', 'blue', 'green', 'orange']

fig = plt.figure(figsize=(20, 12))

# 3D plot of surfaces showing helical structure
ax1 = fig.add_subplot(221, projection='3d')

for i, s in enumerate(s_values):
    filename = f'QI_stellarator_surface_{i+1}_s{s}.dat'
    data = np.loadtxt(filename, comments='#')
    
    x = data[:, 0]
    y = data[:, 1] 
    z = data[:, 2]
    B = data[:, 3]  # Magnetic field strength
    winding = data[:, 4]  # Helical winding parameter
    
    # Plot surface with color based on helical winding to show structure
    scatter = ax1.scatter(x, y, z, c=winding, cmap='viridis', s=2, alpha=0.8, 
                         label=f'Surface s={s}')

# Add field lines to show helical structure
for i in range(min(8, len(field_lines))):
    filename = f'QI_stellarator_fieldline_{i+1}.dat'
    try:
        data = np.loadtxt(filename, comments='#')
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]
        
        # Plot field line with different colors
        color = plt.cm.tab10(i % 10)
        ax1.plot(x, y, z, color=color, linewidth=1, alpha=0.7)
    except:
        pass

# Add coordinate axes
ax1.plot([0, 2.0], [0, 0], [0, 0], 'r-', linewidth=3, label='X-axis')
ax1.plot([0, 0], [0, 2.0], [0, 0], 'g-', linewidth=3, label='Y-axis')
ax1.plot([0, 0], [0, 0], [0, 0.8], 'b-', linewidth=3, label='Z-axis')

ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_zlabel('Z [m]')
ax1.set_title('Quasi-Isodynamic Stellarator\n3D Helical Structure')
ax1.legend()

# 2D cross-section showing helical structure at different angles
ax2 = fig.add_subplot(222)

# Plot cross-sections at different toroidal angles
phi_angles = [0, np.pi/6, np.pi/3, np.pi/2, 2*np.pi/3, 5*np.pi/6, np.pi]
colors_phi = plt.cm.rainbow(np.linspace(0, 1, len(phi_angles)))

for phi, color in zip(phi_angles, colors_phi):
    x_cross = []
    y_cross = []
    
    for s in s_values:
        # Generate cross-section points with true stellarator geometry
        theta = np.linspace(0, 2*np.pi, 200)
        r_base = s * 0.2
        # Strong helical modulation
        helical_mod = 0.25 * np.sin(5 * phi) * 0.3
        r = r_base * (1 + helical_mod)
        R_shift = 0.25 * np.cos(5 * phi) * 0.2
        Z_shift = 0.25 * np.sin(5 * phi) * 0.1
        
        R = 1.0 + r * np.cos(theta) + R_shift
        Z = r * np.sin(theta) + Z_shift
        
        x_cross.extend(R * np.cos(phi))
        y_cross.extend(R * np.sin(phi))
    
    ax2.scatter(x_cross, y_cross, s=1, alpha=0.6, color=color,
               label=f'φ = {phi*180/np.pi:.0f}°')

ax2.set_xlabel('X [m]')
ax2.set_ylabel('Y [m]')
ax2.set_title('Stellarator Cross-Sections\nShowing Helical Winding')
ax2.set_aspect('equal')
ax2.grid(True, alpha=0.3)
ax2.legend()

# 3D view from different angles
ax3 = fig.add_subplot(223, projection='3d')

# Show only the outermost surface for clarity
i = 3  # s = 0.8
filename = f'QI_stellarator_surface_{i+1}_s{s_values[i]}.dat'
data = np.loadtxt(filename, comments='#')
x = data[:, 0]
y = data[:, 1] 
z = data[:, 2]
winding = data[:, 4]

# Plot surface with helical winding color
scatter = ax3.scatter(x, y, z, c=winding, cmap='plasma', s=3, alpha=0.8)

# Add some field lines
for i in range(0, min(4, len(field_lines)), 2):
    filename = f'QI_stellarator_fieldline_{i+1}.dat'
    try:
        data = np.loadtxt(filename, comments='#')
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]
        ax3.plot(x, y, z, 'k-', linewidth=2, alpha=0.8)
    except:
        pass

ax3.set_xlabel('X [m]')
ax3.set_ylabel('Y [m]')
ax3.set_zlabel('Z [m]')
ax3.set_title('Helical Structure Detail\nOuter Surface (s=0.8)')

# Magnetic field strength profile
ax4 = fig.add_subplot(224)

# Plot B field strength along a field line
if len(field_lines) > 0:
    filename = 'QI_stellarator_fieldline_1.dat'
    try:
        data = np.loadtxt(filename, comments='#')
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]
        
        # Calculate B field along the field line
        B_along_line = []
        for j in range(len(x)):
            R = np.sqrt(x[j]**2 + y[j]**2)
            phi = np.arctan2(y[j], x[j])
            Z = z[j]
            
            # Simplified B field calculation
            B_mag = 1.0 + 0.25 * np.cos(5 * phi) * 0.1
            B_along_line.append(B_mag)
        
        ax4.plot(range(len(B_along_line)), B_along_line, 'b-', linewidth=2)
        ax4.set_xlabel('Field Line Position')
        ax4.set_ylabel('|B| [T]')
        ax4.set_title('Magnetic Field Strength\nAlong Field Line')
        ax4.grid(True, alpha=0.3)
    except:
        pass

plt.tight_layout()
plt.savefig('true_QI_stellarator_3d.png', dpi=300, bbox_inches='tight')
plt.savefig('true_QI_stellarator_3d.pdf', bbox_inches='tight')
plt.show()

print("True Quasi-Isodynamic Stellarator visualization saved!")
