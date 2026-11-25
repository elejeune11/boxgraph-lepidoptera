import math
import numpy as np
from scipy import optimize


def _fs_v(x):
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        fs = math.pi * inv * (8.0 / 3.0 + (8.0 / 5.0) * inv2)
        v = math.pi * inv2 * (16.0 / 15.0 + (32.0 / 105.0) * inv2)
        return fs, v
    x2 = x * x
    x3 = x * x2
    at = math.atan(1.0 / x)
    fs = (2.0 * math.pi / x2) * (x + x3 - (x2 * x2 - 1.0) * at)
    v = (math.pi / (2.0 * x3)) * (x + (2.0 / 3.0) * x3 + x3 * x2 - (x2 - 1.0) * ((x2 + 1.0) ** 2) * at)
    return fs, v


def _a_spindle(x):
    if x > 1e3:
        inv = 1.0 / x
        inv2 = inv * inv
        return inv * (8.0 / 3.0 + (8.0 / 15.0) * inv2 - (8.0 / 105.0) * inv2 * inv2)
    eta0 = 2.0 * math.atan(1.0 / x)
    s0 = math.sin(eta0)
    return 2.0 * (eta0 - s0 * math.cos(eta0)) / (s0 * s0)


def run(inputs, settings):
    A_target = float(inputs["A"])
    Kg = float(inputs["K_over_gamma"])
    omega = float(inputs["omega"])
    n_xi = int(settings.get("n_xi", 48))
    n_eta = int(settings.get("n_eta", 24))
    xtol = float(settings.get("xtol", 1e-6))

    u, wu = np.polynomial.legendre.leggauss(n_xi)
    s = 0.5 * (u + 1.0)
    xi = 0.5 * math.pi * (s * s)
    wxi = math.pi * s * 0.5 * wu
    sxi = np.sin(xi)
    cxi = np.cos(xi)
    sxi2 = sxi * sxi
    cxi2 = cxi * cxi
    vn, wv = np.polynomial.legendre.leggauss(n_eta)
    s_eta_u = 0.5 * (vn + 1.0)
    w_eta_u = 0.5 * wv

    x_hom_est = max(1.05, 1.0 + omega if omega < 1.0 else 2.0 * math.sqrt(omega))
    x_grid = [1.01, 1.08, 1.25, 1.6, min(x_hom_est, 12.0), min(1.3 * x_hom_est, 15.0)]
    t_grid = [0.0, 0.01, 0.1, 0.35, 0.7, 0.95, 1.0]

    def solve_for_V(V):
        c = Kg / (V ** (1.0 / 3.0))

        def eval_F(x, t):
            if x < 1.0 or t < 0.0 or t > 1.0:
                return 1e18
            fs, v = _fs_v(x)
            eta0 = 2.0 * math.atan(1.0 / x)
            s0 = math.sin(eta0)
            c0 = math.cos(eta0)
            if t < 1.0:
                Z0 = 1.0 + sxi * c0
                t2N0 = (t * sxi * c0 + 0.5 * (1.0 - t) * Z0) ** 2 + t * sxi2 * (s0 * s0)
                fan = 4.0 * math.pi * np.sum(wxi * (sxi * s0 / (Z0 * Z0)) * (((1.0 - t) ** 2) * cxi2 * (s0 * s0)) / (4.0 * t2N0))
            else:
                fan = 0.0
            if t > 0.0:
                eta = eta0 * s_eta_u
                weta = eta0 * w_eta_u
                seta = np.sin(eta)[None, :]
                ceta = np.cos(eta)[None, :]
                Z = 1.0 + sxi[:, None] * ceta
                t2N = (t * sxi[:, None] * ceta + 0.5 * (1.0 - t) * Z) ** 2 + t * sxi2[:, None] * (seta * seta)
                h_prod = (sxi2[:, None] * seta) / (Z * Z * Z)
                num = 4.0 * cxi2[:, None] + sxi2[:, None] * (seta * seta)
                W = wxi[:, None] * weta[None, :]
                fel = 2.0 * math.pi * np.sum(W * h_prod * (t * t) * num / t2N)
            else:
                fel = 0.0
            return float(c * fel / (v ** (1.0 / 3.0)) + (fs + omega * fan) / (v ** (2.0 / 3.0)))

        best = (1e18, 1.1, 0.5)
        for xg in x_grid:
            for tg in t_grid:
                val = eval_F(xg, tg)
                if val < best[0]:
                    best = (val, xg, tg)
        res = optimize.minimize(
            lambda p: eval_F(p[0], p[1]),
            [best[1], best[2]],
            method="L-BFGS-B",
            bounds=[(1.0, max(20.0, 4.0 * x_hom_est)), (0.0, 1.0)],
            options={"ftol": xtol * 1e-4, "gtol": xtol * 1e-2, "eps": 1e-7},
        )
        res2 = optimize.minimize(
            lambda p: eval_F(p[0], p[1]),
            res.x,
            method="Nelder-Mead",
            options={"xatol": xtol, "fatol": xtol * 1e-3, "maxiter": 300},
        )
        x_opt, t_opt = res2.x if res2.fun < res.fun else res.x
        x_opt = float(max(1.0, x_opt))
        t_opt = float(min(1.0, max(0.0, t_opt)))
        y_opt = float(1.0 / math.sqrt(t_opt)) if t_opt > 1e-16 else 1e8
        _, v_opt = _fs_v(x_opt)
        R_eq = (V / v_opt) ** (1.0 / 3.0)
        A_sp = R_eq * R_eq * _a_spindle(x_opt)
        return x_opt, y_opt, A_sp

    V_sph = (4.0 / 3.0) * math.pi * ((A_target / math.pi) ** 1.5)

    def residual(logV):
        V = math.exp(logV)
        _, _, A_sp = solve_for_V(V)
        return math.log(A_sp / A_target)

    logV_star = optimize.brentq(
        residual,
        math.log(0.02 * V_sph),
        math.log(2.0 * V_sph),
        xtol=max(xtol, 1e-7),
    )
    V_eq = math.exp(logV_star)
    x_eq, y_eq, _ = solve_for_V(V_eq)
    return {"x_eq": x_eq, "y_eq": y_eq, "V_eq": V_eq}
