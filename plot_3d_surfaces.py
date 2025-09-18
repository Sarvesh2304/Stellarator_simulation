import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read surface data
surfaces = []
s_values = [0.2, 0.4, 0.6, 0.8]
colors = ['red', 'blue', 'green', 'orange']

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

for i, s in enumerate(s_values):
    filename = f'surface_{i+1}_s{s}.dat'
    data = np.loadtxt(filename, comments='#')
    
    x = data[:, 0]
    y = data[:, 1] 
    z = data[:, 2]
    
    # Plot surface points
    ax.scatter(x, y, z, c=colors[i], s=2, alpha=0.6, label=f'Surface s={s}')

# Add coordinate axes
ax.plot([0, 1.5], [0, 0], [0, 0], 'r-', linewidth=3, label='X-axis')
ax.plot([0, 0], [0, 1.5], [0, 0], 'g-', linewidth=3, label='Y-axis')
ax.plot([0, 0], [0, 0], [0, 0.5], 'b-', linewidth=3, label='Z-axis')

ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.set_title('3D Stellarator Magnetic Surfaces')
ax.legend()

plt.tight_layout()
plt.savefig('stellarator_3d_surfaces.png', dpi=300, bbox_inches='tight')
plt.savefig('stellarator_3d_surfaces.pdf', bbox_inches='tight')
plt.show()

print("3D visualization saved as 'stellarator_3d_surfaces.png' and '.pdf'")
