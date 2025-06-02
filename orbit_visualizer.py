"""
This script visualizes the trajectory from the most recent run of transfer_optimization.py.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import itertools
import matplotlib.patches as patches

# Constants
earth_radius = 6371e3  # meters

# 1. Locate latest run directory
dirs = sorted(glob.glob('runs/run_*'))
if not dirs:
    raise FileNotFoundError("No run directories found in 'runs/'")
latest = dirs[-1]

# 2. Load optimal trajectory CSV
files = sorted(glob.glob(os.path.join(latest, 'trajectory_*_optimal.csv')))
if not files:
    raise FileNotFoundError(f"No optimal CSV in {latest}")
csv_path = files[-1]

data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
# columns: t, x, y, vx, vy, thrust_flag
t, x, y, vx, vy, thrust = data.T

# 3. Plot continuous trajectory in blue
fig, ax = plt.subplots(figsize=(8,8))
# Draw Earth as solid blue circle
earth = patches.Circle((0, 0), earth_radius, edgecolor='none', facecolor='skyblue', label='Earth')
ax.add_patch(earth)

# Plot trajectory
ax.plot(x, y, color='royalblue', linewidth=1.5, label='Coasting')

# 4. Highlight thrust segments in red
thrust_idx = np.where(thrust == 1)[0]
# Group consecutive indices into segments
groups = [list(g) for _, g in itertools.groupby(enumerate(thrust_idx), key=lambda iv: iv[0] - iv[1])]
for group in groups:
    idxs = [i for _, i in group]
    ax.plot(x[idxs], y[idxs], color='orange', linewidth=1.5, label='Thrusting' if group == groups[0] else "")

# 5. Finalize plot
ax.set_facecolor('white')  # for contrast
#ax.grid(color='gray', linestyle='--', linewidth=0.5)
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('Optimal Trajectory')
ax.legend(loc='upper right')
ax.set_aspect('equal', 'box')
plt.show()
