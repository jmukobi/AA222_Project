import os
import csv
import math
import datetime
import numpy as np

# —— physics constants ——
MU_EARTH    = 3.986e14      # Earth's gravitational parameter, m^3/s^2
G0          = 9.80665       # standard gravity, m/s^2

# —— spacecraft parameters ——
M0          = 1000.0        # kg, initial mass (dry + prop)
T_THRUST    = 1           # N, constant thrust level
ISP         = 3000.0        # s, specific impulse

# —— mission/economic parameters ——
CP          = 5000.0        # $ per kg propellant
R_REVENUE   = 1000.0       # $ per second in GEO

# —— orbit definitions ——
R_EARTH     = 6371e3        # m, Earth radius
R_PERIGEE   = R_EARTH + 200e3    # initial perigee altitude, m
R_APOGEE    = R_EARTH + 35786e3  # GEO apogee altitude, m
R_GEO       = R_APOGEE           # target circular GEO radius

# integration parameters
dt          = 100           # s, fixed time step for RK4
max_steps   = int(1e7)       # safety limit on steps

# golden-section search parameters
ALPHA_MIN   = 0.1           # min fraction of π for half-angle window
ALPHA_MAX   = 0.4           # max fraction
tol         = .1            # convergence tolerance

# simulation counter
sim_counter = 0


def orbital_elements(r, v):
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    h_vec = np.cross(r, v)
    e_vec = (1/MU_EARTH)*((v_norm**2 - MU_EARTH/r_norm)*r - np.dot(r, v)*v)
    e = np.linalg.norm(e_vec)
    a = 1.0/(2.0/r_norm - v_norm**2/MU_EARTH)
    return a, e


def thrust_acceleration(r, v, m, half_angle_window):
    r_norm = np.linalg.norm(r)
    a_grav = -MU_EARTH * r / r_norm**3
    theta = math.atan2(r[1], r[0])
    apogee_angle = math.pi/2
    delta = abs(theta - apogee_angle)
    if delta > math.pi:
        delta = 2*math.pi - delta
    thrust_on = (delta <= half_angle_window)
    if thrust_on:
        v_hat = v / np.linalg.norm(v)
        a_thrust = (T_THRUST / m) * v_hat
    else:
        a_thrust = np.zeros(2)
    return a_grav + a_thrust, thrust_on


def rk4_step(r, v, m, half_angle_window):
    a1, t1 = thrust_acceleration(r, v, m, half_angle_window)
    k1_r = v
    k1_v = a1
    k1_m = -T_THRUST/(ISP*G0) if t1 else 0.0

    r2 = r + 0.5*dt*k1_r
    v2 = v + 0.5*dt*k1_v
    m2 = m + 0.5*dt*k1_m
    a2, t2 = thrust_acceleration(r2, v2, m2, half_angle_window)
    k2_r = v2
    k2_v = a2
    k2_m = -T_THRUST/(ISP*G0) if t2 else 0.0

    r3 = r + 0.5*dt*k2_r
    v3 = v + 0.5*dt*k2_v
    m3 = m + 0.5*dt*k2_m
    a3, t3 = thrust_acceleration(r3, v3, m3, half_angle_window)
    k3_r = v3
    k3_v = a3
    k3_m = -T_THRUST/(ISP*G0) if t3 else 0.0

    r4 = r + dt*k3_r
    v4 = v + dt*k3_v
    m4 = m + dt*k3_m
    a4, t4 = thrust_acceleration(r4, v4, m4, half_angle_window)
    k4_r = v4
    k4_v = a4
    k4_m = -T_THRUST/(ISP*G0) if t4 else 0.0

    r_new = r + (dt/6.0)*(k1_r + 2*k2_r + 2*k3_r + k4_r)
    v_new = v + (dt/6.0)*(k1_v + 2*k2_v + 2*k3_v + k4_v)
    m_new = m + (dt/6.0)*(k1_m + 2*k2_m + 2*k3_m + k4_m)
    return r_new, v_new, m_new, t1


def simulate(half_angle_window, run_dir, optimal=False):
    global sim_counter
    sim_counter += 1
    iteration = sim_counter

    print(f"Starting simulation {iteration}{' (optimal)' if optimal else ''}: half_angle_window={half_angle_window:.4f} rad")
    csv_path = os.path.join(run_dir, f"trajectory_{iteration}{'_optimal' if optimal else ''}.csv")

    r = np.array([0.0, R_APOGEE])
    v0 = math.sqrt(MU_EARTH*(2.0/R_APOGEE - 1.0/((R_PERIGEE + R_APOGEE)/2.0)))
    v = np.array([v0, 0.0])
    m = M0
    t = 0.0

    log = []
    for step in range(max_steps):
        r, v, m, thrust_flag = rk4_step(r, v, m, half_angle_window)
        t += dt
        log.append([t, r[0], r[1], v[0], v[1], int(thrust_flag)])

        a, e = orbital_elements(r, v)
        rp = a*(1.0 - e)
        if step % 100000 == 0:
            print(f"  step {step}, time {t:.1f}s, perigee {rp - R_EARTH:.1f} m above surface")
        if rp >= R_GEO:
            print(f"  reached GEO perigee at step {step}, time {t:.1f}s")
            break

    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["t","x","y","vx","vy","thrust"])
        writer.writerows(log)

    mass_used = M0 - m
    print(f"Simulation {iteration} complete: propellant used={mass_used:.2f} kg, time={t:.1f}s\n")
    return mass_used, t


def cost_function(half_angle_window, run_dir):
    m, t = simulate(half_angle_window, run_dir)
    return CP * m + R_REVENUE * t


def golden_section_search(func, a, b, tol=1e-3, run_dir=None):
    invphi = (math.sqrt(5) - 1) / 2
    invphi2 = (3 - math.sqrt(5)) / 2
    c = a + invphi2*(b - a)
    d = a + invphi*(b - a)
    fc = func(c, run_dir)
    fd = func(d, run_dir)
    iter_count = 1
    print("Starting golden-section search")
    print(f"Iter {iter_count}: a={a:.4f}, b={b:.4f}, c={c:.4f}, d={d:.4f}, fc={fc:.2f}, fd={fd:.2f}")
    while abs(b - a) > tol:
        iter_count += 1
        if fc < fd:
            b, d, fd = d, c, fc
            c = a + invphi2*(b - a)
            fc = func(c, run_dir)
        else:
            a, c, fc = c, d, fd
            d = a + invphi*(b - a)
            fd = func(d, run_dir)
        print(f"Iter {iter_count}: a={a:.4f}, b={b:.4f}, c={c:.4f}, d={d:.4f}, fc={fc:.2f}, fd={fd:.2f}")
    x_opt = (a + b) / 2
    print(f"Golden-section complete in {iter_count} iterations; x_opt={x_opt:.4f}")
    return x_opt, func(x_opt, run_dir)

if __name__ == "__main__":
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = f"runs/run_{now}"
    os.makedirs(run_dir, exist_ok=True)

    def J(phi_frac, run_dir):
        half_ang = phi_frac * math.pi
        return cost_function(half_ang, run_dir)

    phi_opt, J_opt = golden_section_search(J, ALPHA_MIN, ALPHA_MAX, tol, run_dir)
    print(f"Optimal phi_frac: {phi_opt:.4f}, cost: {J_opt:.2f}")

    simulate(phi_opt * math.pi, run_dir, optimal=True)
