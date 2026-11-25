from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def run(inputs, settings):
    # Numerically solves printed SI Eq. (2): d[F-actin]/dt = k_grow * N_fil * (1 - [F-actin]/[T-actin])
    # for [T-actin] satisfying [F-actin](0) = F_0 and [F-actin](t_end) = F_end, without unit conversion.
    k_grow = float(inputs["k_grow"])
    N_fil = float(inputs["N_fil"])
    F_0 = float(inputs["F_0"])
    F_end = float(inputs["F_end"])
    t_end = float(inputs["t_end"])
    dt = float(settings.get("dt", 1.0))

    rate = k_grow * N_fil
    if F_0 + rate * t_end <= F_end or F_end <= F_0:
        return {"T_actin": float("nan")}

    def residual(T):
        sol = solve_ivp(
            lambda t, y: [rate * (1.0 - y[0] / T)],
            (0.0, t_end),
            [F_0],
            method="RK45",
            max_step=dt,
            rtol=1e-8,
            atol=1e-10,
        )
        return float(sol.y[0, -1]) - F_end

    lo = F_end * (1.0 + 1e-9)
    hi = F_end * 1.1
    while residual(hi) < 0.0 and hi < 1e12:
        hi *= 2.0
    if residual(lo) > 0.0 or residual(hi) < 0.0:
        return {"T_actin": float("nan")}

    T_root = brentq(residual, lo, hi, xtol=1e-8, rtol=1e-8)
    return {"T_actin": float(T_root)}
