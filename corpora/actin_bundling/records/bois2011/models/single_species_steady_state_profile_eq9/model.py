import math
import numpy as np
from scipy.optimize import root


def _rk4_half_period_trace(u_peak, s0_hat, Pe, T_half, n_grid):
    """Integrate Eq. (9) from peak (X=0, u=u_peak, u'=0) to valley (X=T_half) using fixed-step RK4."""
    h = T_half / float(n_grid)

    def f_ode(u, up):
        eu = math.exp(min(max(u, -60.0), 60.0))
        upp = u + Pe * (s0_hat - eu / (1.0 + eu))
        return up, upp, eu

    u = float(u_peak)
    up = 0.0
    integ = 0.0
    us = np.empty(n_grid + 1, dtype=float)
    ups = np.empty(n_grid + 1, dtype=float)
    us[0] = u
    ups[0] = up
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

        if 0 < i < n_grid - 1 and prev_up * up_next < 0.0:
            sign_changes += 1
        prev_up = up_next
        u, up = u_next, up_next
        us[i + 1] = u
        ups[i + 1] = up

    xs = np.linspace(0.0, T_half, n_grid + 1)
    return xs, us, ups, integ / T_half, sign_changes


def run(inputs, settings):
    x_eval = float(inputs["x_over_l"])
    Pe = float(inputs["Pe"])
    L_over_l = float(inputs["L_over_l"])
    c0 = float(inputs["c0"])
    n_peaks = int(round(inputs["n_peaks"]))
    x_peak = float(inputs["x_peak_over_l"])
    n_grid = int(round(settings.get("n_grid", 200)))

    period = L_over_l / float(n_peaks)
    T_half = 0.5 * period
    u0 = math.log(c0)
    s0_homog = c0 / (1.0 + c0)

    def residual(params):
        u_peak, s0_hat = params
        _, _, ups, c_mean, _ = _rk4_half_period_trace(u_peak, s0_hat, Pe, T_half, n_grid)
        return np.array([ups[-1], c_mean - c0], dtype=float)

    candidates = [
        (u0 + 2.2, s0_homog - 0.16),
        (u0 + 1.2, s0_homog - 0.09),
        (u0 + 0.6, s0_homog - 0.03),
        (u0 + 3.0, s0_homog - 0.25),
        (u0 + 0.3, s0_homog - 0.01),
    ]
    solved = None
    for guess in candidates:
        res = root(residual, np.array(guess, dtype=float), method="hybr", tol=1e-13)
        if not res.success:
            continue
        u_peak, s0_hat = float(res.x[0]), float(res.x[1])
        xs, us, ups, c_mean, sc = _rk4_half_period_trace(u_peak, s0_hat, Pe, T_half, n_grid)
        if abs(ups[-1]) < 1e-8 and abs(c_mean - c0) < 1e-8 and (us.max() - us.min()) > 1e-3 and sc == 0:
            solved = (xs, us, ups, s0_hat)
            break
    if solved is None:
        raise RuntimeError(f"Could not solve steady profile for Pe={Pe}, n_peaks={n_peaks}")

    xs, us, ups, s0_hat = solved
    # Fold x_eval relative to x_peak onto [-T_half, T_half] by periodicity:
    # u is even around x_peak; du/dx is odd around x_peak (positive to the left, negative to the right).
    dx = ((x_eval - x_peak + T_half) % period) - T_half
    s_dist = abs(dx)
    sign_dx = 1.0 if dx >= 0.0 else -1.0

    u_val = float(np.interp(s_dist, xs, us))
    up_val = float(sign_dx * np.interp(s_dist, xs, ups))
    c_val = math.exp(u_val)
    return {
        "c": float(c_val),
        "v_over_U": float(up_val / Pe),
        "sigma_hat": float(s0_hat + u_val / Pe),
        "active_stress_hat": float(c_val / (1.0 + c_val)),
    }
