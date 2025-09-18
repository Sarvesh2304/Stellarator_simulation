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
ax1.set_title('Quasi-Isodynamic Stellarator\n3D Magnetic Surfaces & Field Lines')
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
ax2.set_title('Stellarator Cross-Sections\nShowing Helical Structure')
ax2.set_aspect('equal')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('quasi_isodynamic_stellarator_3d.png', dpi=300, bbox_inches='tight')
plt.savefig('quasi_isodynamic_stellarator_3d.pdf', bbox_inches='tight')
plt.show()

print("Quasi-isodynamic stellarator visualization saved!")
