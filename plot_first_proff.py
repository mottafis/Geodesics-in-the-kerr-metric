import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

def polar_to_cartesian(df):
    r = df['r']
    phi = df['phi']
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return x, y

files = {
    "direct": 'direct_orbit.txt',
    "retrograde": 'retrograde_orbit.txt'
}

colors = {
    "direct": 'blue',
    "retrograde": 'red'
}

M = 1.0
a = 0.2
E = 0.97
r0 = 10.0
horizon_radius = M + np.sqrt(M**2 - a**2)

def compute_geodesic(a, E, r0):
    if a == 0:
        r_values = np.linspace(r0, horizon_radius * 1.001, 1000)
        phi_values = np.zeros_like(r_values)
    else:
        A = a**2
        B = -2 * M
        C = -(E**2 - 1)
        discriminant = np.sqrt(B**2 - 4 * A * C)
        u_plus = (-B + discriminant) / (2 * A)
        u_minus = (-B - discriminant) / (2 * A)

        def xi_plus(u): return 1 / (u - u_plus + 1e-10)
        def xi_minus(u): return 1 / (u - u_minus + 1e-10)

        def phi_func(r):
            u = 1 / r
            sqrt_arg1 = E**2 * xi_plus(u)**2 + 2 * (M - a**2 * u_plus) * xi_plus(u) - a**2
            sqrt_arg1 = np.clip(sqrt_arg1, 1e-10, None)
            term1 = np.log(np.abs(2 * E * np.sqrt(sqrt_arg1) + 2 * E**2 * xi_plus(u) + 2 * (M - a**2 * u_plus)))

            sqrt_arg2 = E**2 * xi_minus(u)**2 + 2 * (M - a**2 * u_minus) * xi_minus(u) - a**2
            sqrt_arg2 = np.clip(sqrt_arg2, 1e-10, None)
            term2 = np.log(np.abs(2 * E * np.sqrt(sqrt_arg2) + 2 * E**2 * xi_minus(u) + 2 * (M - a**2 * u_minus)))

            return (term1 - term2) / (a * (u_plus - u_minus))

        r_values = np.linspace(r0, horizon_radius * 1.001, 1000)
        phi_values = np.array([phi_func(r) for r in r_values])
        phi_values -= phi_values[0]

    return r_values, phi_values

r_anal, phi_anal = compute_geodesic(a, E, r0)
x_anal = r_anal * np.cos(phi_anal)
y_anal = r_anal * np.sin(phi_anal)

fig, ax = plt.subplots(figsize=(8, 8))

for orbit_type, filepath in files.items():
    df = pd.read_csv(filepath, delim_whitespace=True, comment='#')
    x_num, y_num = polar_to_cartesian(df)
    ax.plot(x_num, y_num, color=colors[orbit_type],
            label=f'{orbit_type.capitalize()} Numerical',
            alpha=0.7, linewidth=1.5)

ax.plot(x_anal, y_anal, 'g--', linewidth=2, label='Analytical (Direct)')

ax.add_patch(patches.Circle((0, 0), horizon_radius,
            edgecolor='black', facecolor='gray', linestyle='--',
            linewidth=1.5, label='Event Horizon'))

ax.set_aspect('equal')
ax.set_xlim(-12, 12)
ax.set_ylim(-12, 12)
ax.set_xlabel('x [M]')
ax.set_ylabel('y [M]')
ax.set_title(f'Geodesics Comparison (a = {a}M, E = {E})')
ax.legend(loc='upper right')
ax.grid(alpha=0.3)

plt.savefig('geodesics_comparison.pdf', bbox_inches='tight')
plt.show()