import math
from scipy.optimize import minimize_scalar


def _f_and_g(a):
    if a < 1e-2:
        # High-precision Taylor series avoiding catastrophic cancellation for small alpha
        a2 = a * a
        f = (a**3) * (1.0 / 3.0 - a2 / 30.0 + (a2 * a2) / 840.0 - (a2**3) / 45360.0)
        g = (a**5) * (2.0 / 15.0 - 11.0 * a2 / 315.0 + 17.0 * (a2 * a2) / 3780.0 - 461.0 * (a2**3) / 1247400.0)
        return f, g
    s = math.sin(a)
    f = s - a * math.cos(a)
    g = f - (1.0 / 3.0) * (s**3)
    return f, g


def _f_hat(a, w):
    f, g = _f_and_g(a)
    return (4.0 * math.pi * f + 2.0 * math.pi * w * g) / ((2.0 * math.pi * g) ** (2.0 / 3.0))


def run(inputs, settings):
    w = float(inputs["omega"])
    n_bracket = int(settings.get("n_bracket", 200))
    xtol = float(settings.get("xtol", 1e-8))

    # Adaptive lower bound so large-omega asymptotic limits (alpha_eq ~ omega^(-1/2)) are always bracketed
    lo_global = min(1e-4, 0.1 / math.sqrt(max(w, 1.0)))
    hi_global = math.pi - 1e-4

    n_log = max(10, n_bracket // 2)
    n_lin = max(10, n_bracket - n_log)
    log_pts = [
        math.exp(math.log(lo_global) + i * (math.log(0.5) - math.log(lo_global)) / (n_log - 1))
        for i in range(n_log)
    ]
    lin_pts = [0.5 + i * (hi_global - 0.5) / (n_lin - 1) for i in range(1, n_lin)]
    grid = log_pts + lin_pts
    vals = [_f_hat(a, w) for a in grid]
    idx = min(range(len(vals)), key=lambda i: vals[i])
    left = grid[max(0, idx - 1)]
    right = grid[min(len(grid) - 1, idx + 1)]

    # Scale xatol relative to bracket width so tiny alpha_eq is resolved to relative precision xtol
    res = minimize_scalar(
        lambda a: _f_hat(a, w),
        bounds=(left, right),
        method="bounded",
        options={"xatol": xtol * max(1e-12, 0.5 * (left + right))},
    )
    a_eq = float(res.x)
    # sin(a) / (1 - cos(a)) = cot(a/2), stable for arbitrarily small a
    aspect = 1.0 / math.tan(0.5 * a_eq)
    return {
        "alpha_eq": a_eq,
        "aspect_ratio": aspect,
        "F_hat_min": float(res.fun),
    }
