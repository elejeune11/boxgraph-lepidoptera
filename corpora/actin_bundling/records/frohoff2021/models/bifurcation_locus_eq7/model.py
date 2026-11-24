import importlib.util
import math
from pathlib import Path
import numpy as np
from scipy.optimize import root

_BT_PATH = Path(__file__).resolve().parent.parent / "bifurcation_threshold_finder_eq7" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_bg_bifurcation_threshold_finder_eq7", _BT_PATH)
_BT_MOD = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BT_MOD)
_make_F = _BT_MOD._make_F


def run(inputs, settings):
    a_start = inputs["a_start"]
    a_end = inputs["a_end"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    Q = inputs["Q"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    ell = inputs["ell"]
    p1_bar = inputs["phi1_bar"]
    p2_bar = inputs["phi2_bar"]
    n_per = max(1, int(round(inputs["n_period"])))
    amp1 = inputs["init_amp_1"]
    amp2 = inputs["init_amp_2"]

    M = int(settings.get("n_modes", 20))
    ds_max = float(settings.get("ds_max", 0.02))
    tol = float(settings.get("tol", 1e-11))
    max_steps = int(settings.get("max_steps", 250))

    F = _make_F(a_Delta, kappa, rho, alpha, ell, p1_bar, p2_bar, n_per, M)
    ratio = (1.0 / Q) * ((rho + alpha) / (rho - alpha))

    def dp_val(u):
        c1, c2 = u[:M], u[M:]
        return 0.5 * float(np.sum(c1 * c1) + ratio * np.sum(c2 * c2))

    u0 = np.zeros(2 * M)
    u0[0] = amp1
    u0[M] = amp2
    s0 = root(lambda u: F(u, a_start), u0, method="hybr", tol=tol)
    if not s0.success or math.sqrt(0.5 * float(np.sum(s0.x**2))) < 1e-4:
        return {"a_crit": None, "norm_at_crit": None}
    u_prev = s0.x

    step_dir = 1.0 if a_end >= a_start else -1.0
    da0 = step_dir * min(abs(a_end - a_start) * 0.05, 0.5 * ds_max)
    s1 = root(lambda u: F(u, a_start + da0), u_prev, method="hybr", tol=tol)
    if not s1.success:
        return {"a_crit": None, "norm_at_crit": None}

    z_prev = np.append(u_prev, a_start)
    z_cur = np.append(s1.x, a_start + da0)
    dp_prev = dp_val(z_prev[:-1])
    dp_cur = dp_val(z_cur[:-1])
    ds = ds_max

    for _ in range(max_steps):
        if dp_prev * dp_cur <= 0.0:
            # Refine exact drift-pitchfork zero along the chord [z_prev, z_cur]
            t_vec = z_cur - z_prev
            t_norm = float(np.linalg.norm(t_vec))
            t_hat = t_vec / t_norm
            lo, hi = 0.0, t_norm
            f_lo = dp_prev
            z_best = z_cur
            for _ in range(45):
                mid = 0.5 * (lo + hi)
                z_pred = z_prev + mid * t_hat
                sol = root(
                    lambda z: np.append(F(z[:-1], z[-1]), float(np.dot(t_hat, z - z_pred))),
                    z_pred,
                    method="hybr",
                    tol=tol,
                )
                if not sol.success:
                    break
                f_mid = dp_val(sol.x[:-1])
                z_best = sol.x
                if abs(f_mid) < 1e-11:
                    break
                if f_lo * f_mid <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    f_lo = f_mid
            return {
                "a_crit": float(z_best[-1]),
                "norm_at_crit": math.sqrt(0.5 * float(np.sum(z_best[:-1] ** 2))),
            }

        t_vec = z_cur - z_prev
        t_hat = t_vec / max(float(np.linalg.norm(t_vec)), 1e-15)
        ok = False
        while ds >= 1e-4:
            z_pred = z_cur + ds * t_hat
            sol = root(
                lambda z: np.append(F(z[:-1], z[-1]), float(np.dot(t_hat, z - z_pred))),
                z_pred,
                method="hybr",
                tol=tol,
            )
            if sol.success and np.max(np.abs(sol.fun)) < 1e-7 and math.sqrt(0.5 * float(np.sum(sol.x[:-1] ** 2))) > 1e-4:
                z_prev, z_cur = z_cur, sol.x
                dp_prev, dp_cur = dp_cur, dp_val(z_cur[:-1])
                ds = min(ds_max, ds * 1.25)
                ok = True
                break
            ds *= 0.5
        if not ok:
            break

    return {"a_crit": None, "norm_at_crit": None}
