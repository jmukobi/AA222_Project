import numpy as np
from scipy.integrate import solve_ivp

# ─── MODULE CONSTANTS ─────────────────────────────────────────────────────────
MU      = 3.986e14        # Earth's gravitational parameter, m^3/s^2
R_P     = 6578e3          # Initial perigee radius (LEO), m
R_A     = 42164e3         # Target apogee radius (GEO), m
T_MAX   = 0.2             # Maximum thrust, N
I_SP    = 1600            # Specific impulse, s
G0      = 9.80665         # Standard gravity, m/s^2
M0      = 500             # Initial spacecraft mass, kg

# Cost weights
C_P     = 1000            # Cost per kg propellant ($/kg)
R_T     = 100             # Cost per second of transfer time ($/s)

# Numeric settings
T_MAX_FACTOR = 100        # factor × orbital period to bound max integration time
NUM_POINTS   = 5000       # points for coast propagation

# ─── ORBIT INITIALIZATION ─────────────────────────────────────────────────────
def init_orbit_parameters():
    a_initial = 0.5 * (R_P + R_A)
    T_est = 2 * np.pi * np.sqrt(a_initial**3 / MU)
    v_perigee = np.sqrt(MU * (2 / R_P - 1 / a_initial))
    y0 = np.array([R_P, 0.0, 0.0, v_perigee, M0])
    return a_initial, T_est, y0

# ─── EVENT FINDERS ─────────────────────────────────────────────────────────────
def find_orbit_events(y0, T_est):
    def coast_odes(t, y):
        x, y_, vx, vy, _ = y
        r = np.hypot(x, y_)
        ax = -MU * x / r**3
        ay = -MU * y_ / r**3
        return [vx, vy, ax, ay, 0.0]

    t_eval = np.linspace(0, T_est, NUM_POINTS)
    sol = solve_ivp(coast_odes, (0, T_est), y0, t_eval=t_eval, max_step=T_est/NUM_POINTS)

    r_hist = np.hypot(sol.y[0], sol.y[1])
    i_peri = np.argmin(r_hist)
    t_peri = sol.t[i_peri]

    vx_pe, vy_pe = sol.y[2, i_peri], sol.y[3, i_peri]
    E = 0.5 * (vx_pe**2 + vy_pe**2) - MU / r_hist[i_peri]
    a_true = -MU / (2 * E)
    T_true = 2 * np.pi * np.sqrt(a_true**3 / MU)
    t_apo = t_peri + T_true / 2
    return t_peri, t_apo, T_true

# ─── PROPAGATION ROUTINES ─────────────────────────────────────────────────────
def propagate_orbit(y0, thrust_schedule, alpha, r_target, t_max):
    def full_odes(t, y):
        x, y_, vx, vy, m = y
        r = np.hypot(x, y_)
        ax = -MU * x / r**3
        ay = -MU * y_ / r**3
        thrust_on = any(start <= t <= end for start, end in thrust_schedule)
        if thrust_on and m > 0:
            ux, uy = x/r, y_/r
            a_t = alpha * T_MAX / m
            ax += a_t * ux
            ay += a_t * uy
            dm_dt = -alpha * T_MAX / (I_SP * G0)
        else:
            dm_dt = 0.0
        return [vx, vy, ax, ay, dm_dt]

    def reach_radius(t, y):
        x, y_, *_ = y
        return np.hypot(x, y_) - r_target
    reach_radius.terminal = True
    reach_radius.direction = 1

    sol = solve_ivp(
        full_odes,
        (0, t_max),
        y0,
        events=reach_radius,
        max_step=t_max/NUM_POINTS
    )
    return sol

# ─── COST FUNCTION ─────────────────────────────────────────────────────────────
def cost_for_alpha(alpha, burn_fraction):
    """
    Compute J = C_P*m_prop + R_T*t_final for throttle fraction alpha
    and burn duration as a fraction of the orbital period.
    """
    a_initial, T_est, y0 = init_orbit_parameters()
    t_peri, t_apo, T_true = find_orbit_events(y0, T_est)
    burn_duration = burn_fraction * T_true
    thrust_schedule = [(t_apo, t_apo + burn_duration)]
    sol = propagate_orbit(y0, thrust_schedule, alpha, R_A, T_true * T_MAX_FACTOR)
    if not sol.t_events[0].size:
        return np.inf
    t_final = sol.t_events[0][0]
    m_final = sol.y_events[0][0][4]
    m_prop = M0 - m_final
    return C_P * m_prop + R_T * t_final

# ─── OPTIMIZATION ROUTINE ─────────────────────────────────────────────────────
def optimize_alpha(burn_fraction):
    # placeholder: return fixed alpha and cost
    class Result:
        pass
    res = Result()
    res.x = 0.5
    res.fun = cost_for_alpha(res.x, burn_fraction)
    return res

# ─── MAIN EXECUTION ───────────────────────────────────────────────────────────
def main():
    burn_fraction = 0.5
    result = optimize_alpha(burn_fraction)
    print(f"Optimal alpha: {result.x:.3f}, cost: {result.fun:.2f}")

if __name__ == '__main__':
    main()
