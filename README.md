# AA222\_PROJECT

This repository contains scripts to optimize a spacecraft transfer to GEO (Geostationary Earth Orbit) using a simple thrust-on/window method and to visualize the resulting trajectory.

## File Structure

```
AA222_PROJECT/
├── runs/                             # Auto-generated run directories storing CSV outputs
│   ├── run_YYYYMMDD_HHMMSS/          # Example: run_20250516_174343
│   │   ├── trajectory_1.csv
│   │   ├── trajectory_2.csv
│   │   ├── ...
│   │   └── trajectory_7_optimal.csv  # Optimal trajectory from golden-section search
│   └── ...
├── transfer_optimization.py         # Main optimization script
└── orbit_visualizer.py              # Script to plot the optimal trajectory
```

## Dependencies

* Python 3.x
* NumPy
* Matplotlib

You can install the required Python packages using pip:

```bash
pip install numpy matplotlib
```

## Scripts Overview

### `transfer_optimization.py`

This is the main script that:

1. Defines spacecraft and mission parameters (e.g., thrust level, specific impulse, cost parameters).
2. Simulates an orbital transfer from a 200 km perigee to GEO using a thrust window around apogee.
3. Uses a golden-section search to find the optimal half-angle window fraction (`phi_frac`) that minimizes a cost function:

   * **Cost function** = $C_P \times \text{(propellant used)} + R_{REVENUE} \times \text{(transfer time)}$
4. Logs the spacecraft trajectory (time, position, velocity, thrust flag) into CSV files under `runs/run_YYYYMMDD_HHMMSS/`.
5. Outputs the optimal `phi_frac` and its associated cost, then simulates and saves the final optimal trajectory.

#### Usage

1. Run the script:

   ```bash
   python transfer_optimization.py
   ```

2. The script will create a new run directory under `runs/` named `run_<timestamp>/`.

3. Trajectory CSV files (e.g., `trajectory_1.csv`, ..., `trajectory_7_optimal.csv`) will be saved in the new run folder.

4. Check the console output to see progress details (perigee altitude updates, cost evaluations, optimal window fraction).

#### Key Parameters (Editable in Script)

* **`MU_EARTH`**: Earth's gravitational parameter (m³/s²).
* **`G0`**: Standard gravity (m/s²).
* **`M0`**: Initial spacecraft mass (kg).
* **`T_THRUST`**: Constant thrust level (N).
* **`ISP`**: Specific impulse (s).
* **`CP`**: Cost per kg of propellant (\$/kg).
* **`R_REVENUE`**: Revenue rate (\$/s) while in GEO.
* **`R_PERIGEE`**: Initial perigee radius (m).
* **`R_APOGEE`**: Initial apogee radius (m).
* **`dt`**: Time step for RK4 integration (s).
* **`ALPHA_MIN`**, **`ALPHA_MAX`**: Fraction bounds (of π) for golden-section search window.
* **`tol`**: Convergence tolerance for golden-section search.

Adjust these constants at the top of `transfer_optimization.py` to experiment with different mission/vehicle scenarios.

### `orbit_visualizer.py`

This script loads the latest optimal trajectory CSV from the most recent run directory and generates a 2D plot:

* Draws Earth as a circle with radius \~6371 km.
* Plots the spacecraft trajectory in blue (coasting) and highlights thrust segments in orange.

#### Usage

1. Ensure you have at least one completed run directory under `runs/`.

2. Run the visualizer:

   ```bash
   python orbit_visualizer.py
   ```

3. A plot window will appear showing the trajectory in an Earth-centered frame.

## Example Workflow

1. Run the optimization:

   ```bash
   python transfer_optimization.py
   ```

   * Monitor console output for progress and final optimal window fraction.

2. Plot the optimal trajectory:

   ```bash
   python orbit_visualizer.py
   ```

   * View the plotted orbit and thrust segments.

## CSV Output Format

Each `trajectory_*.csv` file contains rows with the following columns:

```
t, x, y, vx, vy, thrust
```

* **`t`**: Simulation time (s).
* **`x`, `y`**: Spacecraft position in the orbital plane (m).
* **`vx`, `vy`**: Spacecraft velocity components (m/s).
* **`thrust`**: Flag (0 or 1) indicating whether thrust is on.
