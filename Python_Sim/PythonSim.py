import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from numpy.linalg import LinAlgError, solve
import time
start_time = time.perf_counter()

t_max = 900.0 # seconds
NUM_TIME_STEPS = 500
time_step = t_max / NUM_TIME_STEPS
t = np.linspace(0, t_max, NUM_TIME_STEPS)

NUM_NODES_PER_CM = 4
plates_config = [ 
    {"name": "Load (MC)", "temp": 0,  "length_below": 0.0},   
    {"name": "Plate 1",   "temp": 0.010, "length_below": 0.15}, 
    {"name": "Plate 2",   "temp": 0.10,   "length_below": 0.15}, 
    {"name": "Plate 3",   "temp": 1.0,   "length_below": 0.15}, 
    {"name": "Plate 4",   "temp": 4.0,  "length_below": 0.15}, 
    {"name": "Plate 5",   "temp": 50.0,  "length_below": 0.40}, 
    {"name": "Plate 6",   "temp": 300.0, "length_below": 0.40}, 
]

LOAD_RESISTANCE = 1.0 
CURRENT_AMPS = 0.0
WIRE_DIAMETER = 0.001 # wasn't used, but could be to reason cross sectional area, I just made unity.

current_pos = 0
for p in plates_config:
    current_pos += p["length_below"]
    p["abs_position"] = current_pos

LENGTH = plates_config[-1]["abs_position"]
NUM_NODES = int(LENGTH * 100 * NUM_NODES_PER_CM)
dx = LENGTH / NUM_NODES
x_coords = np.linspace(0, LENGTH, NUM_NODES)

MATERIALS = {
    'cu_wire_properties' : {
    'density': 8960,  # kg/m^3
    'specific_heat': 385,  # J/(kg·K)
    'thermal_conductivity': 401,  # W/(m·K)
    'rho_elec': lambda T: 1.68e-8 if T < 20 else 1.68e-8 * (1 + 0.0039 * (T - 20)),
    'diffusion_coeff': 401 / (8960 * 385)  # k / (rho * Cp)
    }
}

T = np.zeros((NUM_NODES, NUM_TIME_STEPS))
T[:, 0] = 300.0  # Initial temperature in Kelvin
fixed_nodes = []
for p in plates_config:
    if p["temp"] is not None:
        idx = np.abs(x_coords - p["abs_position"]).argmin()
        fixed_nodes.append({"idx": idx, "temp": p["temp"], "name": p["name"]})
        T[idx, :] = p["temp"] # Note: Fixed temperature set for all time steps


lam = (MATERIALS['cu_wire_properties']['diffusion_coeff'] * time_step) / (2 * dx**2)
N_inner = NUM_NODES - 2

# Matrix A (LHS: Coefficients for T^(j+1))
main_diag_A = np.full(N_inner, 2 + 2 * lam)
off_diag_A = np.full(N_inner - 1, -lam)
A = np.diag(main_diag_A, 0) + np.diag(off_diag_A, 1) + np.diag(off_diag_A, -1)

# Matrix B (RHS: Coefficients for T^j)
main_diag_B = np.full(N_inner, 2 - 2 * lam)
off_diag_B = np.full(N_inner - 1, lam)
B = np.diag(main_diag_B, 0) + np.diag(off_diag_B, 1) + np.diag(off_diag_B, -1)


for j in range(0, NUM_TIME_STEPS - 1):
    b = T[1: -1, j].copy()
    b = np.dot(B, b)
    
    b[0] += lam * T[0, j]
    b[-1] += lam * T[NUM_NODES - 1, j+1]

    solution = np.linalg.solve(A, b)
    T[1: -1, j + 1] = solution

    for node in fixed_nodes:
        T[node["idx"], j + 1] = node["temp"]

end_time = time.perf_counter()

print(f"Runtime: {end_time - start_time} seconds")

# Plotting _________________

R = np.linspace(1, 0 , NUM_TIME_STEPS)
B = np.linspace(0, 1 , NUM_TIME_STEPS)
G = 0

plt.figure(figsize=(12, 6))
for i in range(NUM_TIME_STEPS//50 -1):
    plt.plot(x_coords, T[:, (i+1)*50], color=(R[i*50], G, B[i*50]), label=f'Time={(i+1)*50*time_step:.1f}s')
    
plt.axhline(y=9.2, color='purple', linestyle='--', label='9.2 K (NbTi Tc)')
plt.xlabel('Position along wire (m)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Distribution along the Wire over Time')

plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    borderaxespad=0.,
    fontsize=8
)

plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()



Q_LEAK_6 = CURRENT_AMPS*CURRENT_AMPS*LOAD_RESISTANCE # Heat leak into Plate 6
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(12, 6))
cmap = plt.get_cmap('turbo')
norm = plt.Normalize(vmin=np.min(T), vmax=np.max(T))

points = np.array([x_coords, np.zeros_like(x_coords)]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap=cmap, norm=norm)
lc.set_array(T[:, -1])
lc.set_linewidth(4)
ax.add_collection(lc)

block_size = LENGTH * 0.04
for p in plates_config:
    x_pos = p["abs_position"]
    local_temp = T[np.abs(x_coords - x_pos).argmin(), -1]
    color = cmap(norm(local_temp))
    
    rect = patches.Rectangle(
        (x_pos - block_size/2, -block_size/2), 
        block_size, block_size,                
        linewidth=1, edgecolor='lime', facecolor=color, zorder=10
    )
    ax.add_patch(rect)
    
    ax.add_patch(patches.Rectangle(
        (x_pos - block_size/2, -block_size/2),
        block_size, block_size,
        linewidth=0, facecolor='none', edgecolor='black', hatch='XX', alpha=0.3, zorder=11
    ))

    label_text = f"{p['name']}\n{local_temp:.2f} K"
    ax.text(x_pos, -block_size*0.9, label_text, color='white', ha='center', va='top', fontsize=8)

ax.text(0, block_size, f"No Load\n(Q=0 W)", color='lime', ha='center', fontweight='bold')
ax.text(LENGTH, block_size, f"Q leak: {Q_LEAK_6} W", color='yellow', ha='center', fontweight='bold')


ax.set_xlim(-block_size, LENGTH + block_size)
ax.set_ylim(-block_size*2.5, block_size*2)
ax.set_aspect('equal')
ax.axis('off')

# Arrow pointing to 9.2K (Critical Temperature of NbTi)
crit_temp = 9.2
crit_idx = np.abs(T[:, -1] - crit_temp).argmin()
ax.annotate('9.2 K (NbTi Tc)', xy=(x_coords[crit_idx], block_size/2), xytext=(x_coords[crit_idx], block_size),
            arrowprops=dict(facecolor='cyan', shrink=0.05),
            bbox=dict(boxstyle="round,pad=0.3", fc="cyan", ec="none", alpha=0.8),
            color='black', fontweight='bold', ha='center')

# Text with Location of 9.2K (meters from load)
ax.text((LENGTH)/2, -block_size*2.5, 
        f"Location of 9.2 K: {x_coords[crit_idx]:.3f} m", 
        color='white', fontsize=10, ha='center', va='center')

ax.text((LENGTH)/2, -block_size*3.0, 
        f"Distance from 4 K Plate: {x_coords[crit_idx] - plates_config[4]['abs_position']:.3f} m", 
        color='white', fontsize=10, ha='center', va='center')


cbar = plt.colorbar(lc, ax=ax, orientation='horizontal', pad=0.2, fraction=0.05)
cbar.set_label('Temperature (K)', color='white')
cbar.ax.xaxis.set_tick_params(color='white')
plt.setp(plt.getp(cbar.ax.axes, 'xticklabels'), color='white')
plt.title(f"Thermal Model: Zero-Load Conduction Leak (I=0A)", color='white', fontsize=14, pad=10)
plt.show()

# ==========================================
# TODO: Implement the NbTi portion
# Determine the maximum current the wire can carry without exceeding 9.2 K
# Determine the resistance of the wire at various temperatures

# Note: implementing NbTi requires additional material properties and logic not present in the original code.
# That logic consists of checking the temperature profile for exceeding 9.2 K and adjusting current accordingly.
# And defining spatially varying resistivity (since some portions of the wire may be NbTi and others Cu).