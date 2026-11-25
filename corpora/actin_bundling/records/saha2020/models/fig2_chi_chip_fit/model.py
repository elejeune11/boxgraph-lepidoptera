import math
import numpy as np
from scipy.optimize import brentq, minimize


def _fn_vec(p1, p2, c11, c12, c21, c22, chi, chip, alpha, is_ci):
    d11 = -(2.0 * (p1 - c11) ** 2 + 8.0 * (p1 - c11) * (p1 - c12) + 2.0 * (p1 - c12) ** 2 + 2.0 * p2 * p2 * chip)
    d22 = -(2.0 * (p2 - c21) ** 2 + 8.0 * (p2 - c21) * (p2 - c22) + 2.0 * (p2 - c22) ** 2 + 2.0 * p1 * p1 * chip)
    d12 = -((chi + alpha) + 4.0 * p1 * p2 * chip)
    d21 = -((chi - alpha) + 4.0 * p1 * p2 * chip)
    disc = (d11 - d22) ** 2 + 4.0 * d12 * d21
    det = d11 * d22 - d12 * d21
    return np.where(is_ci, disc, det)


def _fn_scalar(p1, p2, c11, c12, c21, c22, chi, chip, alpha, kind):
    d11 = -(2.0 * (p1 - c11) ** 2 + 8.0 * (p1 - c11) * (p1 - c12) + 2.0 * (p1 - c12) ** 2 + 2.0 * p2 * p2 * chip)
    d22 = -(2.0 * (p2 - c21) ** 2 + 8.0 * (p2 - c21) * (p2 - c22) + 2.0 * (p2 - c22) ** 2 + 2.0 * p1 * p1 * chip)
    d12 = -((chi + alpha) + 4.0 * p1 * p2 * chip)
    d21 = -((chi - alpha) + 4.0 * p1 * p2 * chip)
    disc = (d11 - d22) ** 2 + 4.0 * d12 * d21
    if kind == "ci":
        return disc
    tr = d11 + d22
    return 0.5 * (tr + math.sqrt(disc)) if disc >= 0.0 else 0.5 * tr


def _root(p1, c11, c12, c21, c22, chi, chip, alpha, kind, branch, y_lo, y_hi):
    ys = np.linspace(y_lo, y_hi, 120)
    vals = [_fn_scalar(p1, float(y), c11, c12, c21, c22, chi, chip, alpha, kind) for y in ys]
    roots = []
    for i in range(len(ys) - 1):
        if vals[i] * vals[i + 1] <= 0.0:
            try:
                r = brentq(lambda y: _fn_scalar(p1, y, c11, c12, c21, c22, chi, chip, alpha, kind), float(ys[i]), float(ys[i + 1]))
                roots.append(float(r))
            except Exception:
                pass
    if not roots:
        return None
    return max(roots) if branch == "upper" else min(roots)


def _eval_curve(obs, c11, c12, c21, c22, chi, chip, alpha, kind, branch, y_lo, y_hi):
    xs = list(obs["x"])
    ys_out = []
    for x in xs:
        r = _root(float(x), c11, c12, c21, c22, chi, chip, alpha, kind, branch, y_lo, y_hi)
        ys_out.append(float("nan") if r is None else r)
    return {"x": xs, "y": ys_out}


def _saddle_alpha_split(c11, c12, c21, c22, chi, chip):
    def det0(xy):
        p1, p2 = float(xy[0]), float(xy[1])
        d11 = -(2.0 * (p1 - c11) ** 2 + 8.0 * (p1 - c11) * (p1 - c12) + 2.0 * (p1 - c12) ** 2 + 2.0 * p2 * p2 * chip)
        d22 = -(2.0 * (p2 - c21) ** 2 + 8.0 * (p2 - c21) * (p2 - c22) + 2.0 * (p2 - c22) ** 2 + 2.0 * p1 * p1 * chip)
        s = chi + 4.0 * p1 * p2 * chip
        return d11 * d22 - s * s

    def grad_norm(xy):
        h = 1e-6
        gx = (det0((xy[0] + h, xy[1])) - det0((xy[0] - h, xy[1]))) / (2 * h)
        gy = (det0((xy[0], xy[1] + h)) - det0((xy[0], xy[1] - h))) / (2 * h)
        return gx * gx + gy * gy

    res = minimize(grad_norm, [0.18, 0.14], method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-13})
    return math.sqrt(max(0.0, -det0(res.x)))


def run(inputs, settings):
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    alpha_a = float(inputs["alpha_a"])
    alpha_b = float(inputs["alpha_b"])

    curves = [
        ("cr_2a_upper", inputs["cr_2a_upper_obs"], alpha_a, "cr", "upper", -0.5, 0.6),
        ("cr_2a_lower", inputs["cr_2a_lower_obs"], alpha_a, "cr", "lower", -0.5, 0.5),
        ("ci_2a_lower", inputs["ci_2a_lower_obs"], alpha_a, "ci", "lower", 0.05, 0.65),
        ("cr_2b_ul", inputs["cr_2b_ul_obs"], alpha_b, "cr", "upper", 0.05, 0.55),
        ("cr_2b_lr", inputs["cr_2b_lr_obs"], alpha_b, "cr", "lower", -0.5, 0.15),
        ("ci_2b_ll", inputs["ci_2b_ll_obs"], alpha_b, "ci", "upper", -0.65, -0.05),
        ("cr_2c_0", inputs["cr_2c_0_obs"], 0.0, "cr", "upper", -0.45, 0.6),
        ("cr_2c_09_up", inputs["cr_2c_09_up_obs"], 0.9, "cr", "upper", 0.2, 0.45),
        ("cr_2c_09_lo", inputs["cr_2c_09_lo_obs"], 0.9, "cr", "lower", 0.2, 0.45),
    ]

    # Stack all points into flat arrays for vectorized Newton-distance optimization (~5 ms total!)
    p1_all, p2_all, al_all, ci_all = [], [], [], []
    for _, obs, al, kind, _, _, _ in curves:
        for x, y in zip(obs["x"], obs["y"]):
            p1_all.append(float(x))
            p2_all.append(float(y))
            al_all.append(float(al))
            ci_all.append(kind == "ci")
    p1_a = np.array(p1_all)
    p2_a = np.array(p2_all)
    al_a = np.array(al_all)
    ci_a = np.array(ci_all, dtype=bool)
    eps = 1e-5

    def loss(params):
        chi, chip = float(params[0]), float(params[1])
        f0 = _fn_vec(p1_a, p2_a, c11, c12, c21, c22, chi, chip, al_a, ci_a)
        df_dy = (_fn_vec(p1_a, p2_a + eps, c11, c12, c21, c22, chi, chip, al_a, ci_a)
                 - _fn_vec(p1_a, p2_a - eps, c11, c12, c21, c22, chi, chip, al_a, ci_a)) / (2.0 * eps)
        df_dx = (_fn_vec(p1_a + eps, p2_a, c11, c12, c21, c22, chi, chip, al_a, ci_a)
                 - _fn_vec(p1_a - eps, p2_a, c11, c12, c21, c22, chi, chip, al_a, ci_a)) / (2.0 * eps)
        grad = np.hypot(df_dx, df_dy) + 1e-4
        dist = f0 / grad
        return float(np.mean(dist * dist))

    res = minimize(loss, [-0.50, 0.85], method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-10})
    chi_fit, chip_fit = float(res.x[0]), float(res.x[1])
    alpha_split = _saddle_alpha_split(c11, c12, c21, c22, chi_fit, chip_fit)

    out = {
        "chi_fit": chi_fit,
        "chip_fit": chip_fit,
        "alpha_split_fit": alpha_split,
    }
    for name, obs, al, kind, br, ylo, yhi in curves:
        out[f"{name}_fit"] = _eval_curve(obs, c11, c12, c21, c22, chi_fit, chip_fit, al, kind, br, ylo, yhi)
    return out
