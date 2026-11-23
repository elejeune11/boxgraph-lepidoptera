import math
import numpy as np
from scipy.optimize import root


def _rk4_half_period(u_peak, s0_hat, Pe, T_half, n_grid):
    """Integrate Eq. (9) from peak (X=0, u=u_peak, u'=0) to valley (X=T_half) using fixed-step RK4."""
    h = T_half / float(n_grid)

    def f_ode(u, up):
        eu = math.exp(min(max(u, -60.0), 60.0))
        upp = u + Pe * (s0_hat - eu / (1.0 + eu))
        return up, upp, eu

    u = float(u_peak)
    up = 0.0
    integ = 0.0
    u_min = u
    u_max = u
    up_abs_max = 0.0
    # Track monotonicity: starting from a peak (u_peak > ln(c0)), u' should be <= 0 on (0, T_half)
    sign_changes = 0
    prev_up = 0.0

    for i in range(n_grid):
        k1_u, k1_up, k1_e = f_ode(u, up)
        k2_u, k2_up, k2_e = f_ode(u + 0.5 * h * k1_u, up + 0.5 * h * k1_up)
        k3_u, k3_up, k3_e = f_ode(u + 0.5 * h * k2_u, up + 0.5 * h * k2_up)
        k4_u, k4_up, k4_e = f_ode(u + h * k3_u, up + h * k3_up)

        u_next = u + (h / 6.0) * (k1_u + 2.0 * k2_u + 2.0 * k3_u + k4_u)
        up_next = up + (h / 6.0) * (k1_up + 2.0 * k2_up + 2.0 * k3_up + k4_up)
        integ += (h / 6.0) * (k1_e + 2.0 * k2_e + 2.0 * k3_e + k4_e)

        if i > 0 and i < n_grid - 1 and prev_up * up_next < 0.0:
            sign_changes += 1
        prev_up = up_next
        u, up = u_next, up_next
        if u < u_min:
            u_min = u
        if u > u_max:
            u_max = u
        if abs(up) > up_abs_max:
            up_abs_max = abs(up)

    return u, up, integ / T_half, u_min, u_max, up_abs_max, sign_changes


def run(inputs, settings):
    Pe = float(inputs["Pe"])
    L_over_l = float(inputs["L_over_l"])
    c0 = float(inputs["c0"])
    n_peaks = int(round(inputs["n_peaks"]))
    n_grid = int(settings.get("n_grid", 100))

    if n_peaks < 1:
        raise ValueError("n_peaks must be >= 1")

    T_half = L_over_l / (2.0 * n_peaks)
    u0 = math.log(c0)
    s0_homog = c0 / (1.0 + c0)

    def residual(params):
        u_peak, s0_hat = params
        _, up_end, c_mean, _, _, _, _ = _rk4_half_period(u_peak, s0_hat, Pe, T_half, n_grid)
        return np.array([up_end, c_mean - c0], dtype=float)

    # Multi-start shooting from peak u_peak > ln(c0) to ensure we land on the fundamental
    # nonhomogeneous branch of period L/n_peaks (no interior zero crossings of u').
    candidates = [
        (u0 + 2.2, s0_homog - 0.16),
        (u0 + 1.2, s0_homog - 0.09),
        (u0 + 0.6, s0_homog - 0.03),
        (u0 + 3.0, s0_homog - 0.25),
        (u0 + 0.3, s0_homog - 0.01),
    ]
    best = None
    for guess in candidates:
        res = root(residual, np.array(guess, dtype=float), method="hybr", tol=1e-13)
        if not res.success:
            continue
        u_peak, s0_hat = float(res.x[0]), float(res.x[1])
        u_end, up_end, c_mean, u_min, u_max, up_abs_max, sc = _rk4_half_period(
            u_peak, s0_hat, Pe, T_half, n_grid
        )
        if abs(up_end) < 1e-8 and abs(c_mean - c0) < 1e-8 and (u_max - u_min) > 1e-3 and sc == 0:
            best = (u_peak, s0_hat, u_min, u_max, c_mean, up_abs_max)
            break

    if best is None:
        raise RuntimeError(
            f"Could not locate nonhomogeneous {n_peaks}-peak steady state for Pe={Pe}, L/l={L_over_l}, c0={c0}"
        )

    _, s0_hat, u_min, u_max, c_mean, up_abs_max = best
    c_max = math.exp(u_max)
    c_min = math.exp(u_min)
    return {
        "c_max": float(c_max),
        "c_min": float(c_min),
        "c_mean": float(c_mean),
        "v_over_U_max": float(up_abs_max / Pe),
        "sigma0_hat": float(s0_hat),
        "sigma_hat_max": float(s0_hat + u_max / Pe),
        "sigma_hat_min": float(s0_hat + u_min / Pe),
        "active_stress_hat_max": float(c_max / (1.0 + c_max)),
    }
