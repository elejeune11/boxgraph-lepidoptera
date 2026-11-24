import math
import numpy as np
from scipy.optimize import brentq, curve_fit, minimize_scalar


def _f_and_g(a):
    if a < 1e-2:
        a2 = a * a
        f = (a**3) * (1.0 / 3.0 - a2 / 30.0 + (a2 * a2) / 840.0 - (a2**3) / 45360.0)
        g = (a**5) * (2.0 / 15.0 - 11.0 * a2 / 315.0 + 17.0 * (a2 * a2) / 3780.0 - 461.0 * (a2**3) / 1247400.0)
        return f, g
    s = math.sin(a)
    f = s - a * math.cos(a)
    g = f - (1.0 / 3.0) * (s**3)
    return f, g


def _fp_and_gp(a):
    if a < 1e-2:
        a2 = a * a
        fp = (a**2) * (1.0 - a2 / 6.0 + (a2 * a2) / 120.0 - (a2**3) / 5040.0)
        gp = (a**4) * (2.0 / 3.0 - 11.0 * a2 / 45.0 + 17.0 * (a2 * a2) / 420.0 - 461.0 * (a2**3) / 113400.0)
        return fp, gp
    s = math.sin(a)
    c = math.cos(a)
    fp = a * s
    gp = a * s - (s**2) * c
    return fp, gp


def _f_hat(a, w):
    f, g = _f_and_g(a)
    return (4.0 * math.pi * f + 2.0 * math.pi * w * g) / ((2.0 * math.pi * g) ** (2.0 / 3.0))


def _find_alpha_eq(w):
    lo_global = min(1e-4, 0.1 / math.sqrt(max(w, 1.0)))
    hi_global = math.pi - 1e-4
    grid = [
        math.exp(math.log(lo_global) + i * (math.log(0.5) - math.log(lo_global)) / 99)
        for i in range(100)
    ]
    grid += [0.5 + i * (hi_global - 0.5) / 99 for i in range(1, 100)]
    vals = [_f_hat(a, w) for a in grid]
    idx = min(range(len(vals)), key=lambda i: vals[i])
    left = grid[max(0, idx - 1)]
    right = grid[min(len(grid) - 1, idx + 1)]
    res = minimize_scalar(
        lambda a: _f_hat(a, w),
        bounds=(left, right),
        method="bounded",
        options={"xatol": 1e-12 * max(1e-12, 0.5 * (left + right))},
    )
    return float(res.x)


def run(inputs, settings):
    w = float(inputs["omega"])
    v_cap = float(inputs["v_cap"])
    L_eq = float(inputs["L_eq"])
    L_i = float(inputs["L_i"])
    t_max = float(inputs["t_max"])
    n_steps = int(settings.get("n_steps", 400))

    a_eq = _find_alpha_eq(w)
    _, g_eq = _f_and_g(a_eq)
    R_eq = L_eq / (2.0 * math.sin(a_eq))
    V = 2.0 * math.pi * (R_eq**3) * g_eq

    def L_of_a(a):
        _, g = _f_and_g(a)
        R = (V / (2.0 * math.pi * g)) ** (1.0 / 3.0)
        return 2.0 * R * math.sin(a)

    # Find initial alpha_i corresponding to L(alpha_i) = L_i at constant V
    if abs(L_i - L_eq) < 1e-12:
        a_i = a_eq
    elif L_i > L_eq:
        lo_a = min(1e-4, 0.25 * a_eq)
        while L_of_a(lo_a) < L_i and lo_a > 1e-10:
            lo_a *= 0.25
        a_i = brentq(lambda a: L_of_a(a) - L_i, lo_a, a_eq, xtol=1e-12)
    else:
        a_i = brentq(lambda a: L_of_a(a) - L_i, a_eq, math.pi - 1e-6, xtol=1e-12)

    two_pi_13 = (2.0 * math.pi) ** (1.0 / 3.0)
    V_23 = V ** (2.0 / 3.0)

    def da_dt(a):
        f, g = _f_and_g(a)
        fp, gp = _fp_and_gp(a)
        dF_da_over_gI = two_pi_13 * V_23 * (
            2.0 * fp * (g ** (-2.0 / 3.0))
            - (4.0 / 3.0) * f * (g ** (-5.0 / 3.0)) * gp
            + (1.0 / 3.0) * w * (g ** (-2.0 / 3.0)) * gp
        )
        # Printed SI p. 9 col 1: dF/dt + D = 0 with D = (1/2) V eta [((1/R) dR/dt)^2 + ((1/a) da/dt)^2]
        mob = 0.5 * V * ((gp / (3.0 * g)) ** 2 + 1.0 / (a**2))
        return -v_cap * dF_da_over_gI / mob

    # Fixed-step RK4 integration so refinements over n_steps directly test time-step convergence
    dt = t_max / float(n_steps)
    ts = np.linspace(0.0, t_max, n_steps + 1)
    alphas = np.empty(n_steps + 1, dtype=float)
    alphas[0] = a_i
    a_curr = a_i
    for i in range(n_steps):
        k1 = da_dt(a_curr)
        k2 = da_dt(a_curr + 0.5 * dt * k1)
        k3 = da_dt(a_curr + 0.5 * dt * k2)
        k4 = da_dt(a_curr + dt * k3)
        a_curr = a_curr + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        alphas[i + 1] = a_curr

    Ls = np.array([L_of_a(float(a)) for a in alphas])

    # Analytical linearized timescale near alpha_eq
    eps = 1e-5 * a_eq
    d2F_hat = (_f_hat(a_eq + eps, w) - 2.0 * _f_hat(a_eq, w) + _f_hat(a_eq - eps, w)) / (eps**2)
    d2F_over_gI = d2F_hat * V_23
    _, gp_eq = _fp_and_gp(a_eq)
    mob_eq = 0.5 * V * ((gp_eq / (3.0 * g_eq)) ** 2 + 1.0 / (a_eq**2))
    tau_lin = float(mob_eq / (v_cap * d2F_over_gI))

    if abs(L_i - L_eq) < 1e-12:
        tau_fit = tau_lin
    else:
        popt, _ = curve_fit(
            lambda t, tau: L_eq + (L_i - L_eq) * np.exp(-t / tau),
            ts,
            Ls,
            p0=[max(1e-3, tau_lin)],
            bounds=(1e-6, 1e8),
        )
        tau_fit = float(popt[0])

    return {
        "tau": tau_fit,
        "tau_linearized": tau_lin,
        "alpha_eq": a_eq,
        "alpha_min": float(np.min(alphas)),
        "L_final": float(Ls[-1]),
    }
