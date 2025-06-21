import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

def polar_to_cartesian(df):
    """Convert (r, phi) to Cartesian coordinates using Kerr-Schild transformation"""
    r = df['r']
    phi = df['phi']
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return x, y

# Configuration
files = {
    "direct": 'direct_orbit.txt',
    "retrograde": 'retrograde_orbit.txt'
}

colors = {
    "direct": 'blue',
    "retrograde": 'red'
}

# Physical parameters
M = 1.0
a = 0.99
horizon_radius = M + np.sqrt(M**2 - a**2)  # Event horizon calculation

# Create plot
fig, ax = plt.subplots(figsize=(7, 7))

# Plot trajectories
for orbit_type, filepath in files.items():
    df = pd.read_csv(filepath, delim_whitespace=True, comment='#')
    x, y = polar_to_cartesian(df)
    ax.plot(x, y, color=colors[orbit_type], 
            label=f'{orbit_type.capitalize()} Orbit', 
            alpha=0.7, linewidth=1.5)

# Add event horizon
horizon = patches.Circle((0, 0), horizon_radius, 
                         edgecolor='black', facecolor='gray',
                         linestyle='--', linewidth=2,
                         label='Event Horizon')
ax.add_patch(horizon)

# Plot aesthetics
ax.set_aspect('equal')
ax.set_xlabel('x [M]')
ax.set_ylabel('y [M]')
ax.set_title('Kerr Black Hole Geodesics (a = 0.99M)')
ax.legend(loc='upper right')
ax.grid(True, alpha=0.3)


plt.savefig('kerr_geodesics.pdf', bbox_inches='tight')
plt.show()
