import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq, minimize_scalar

_EULER_GAMMA = 0.5772156649015329


def _log1pexp(x):
    if x > 50.0:
        return x
    if x < -50.0:
        return math.exp(x)
    return math.log1p(math.exp(x))


def _sigmoid(x):
    if x >= 0.0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def run(inputs, settings):
    c_l_target = inputs["c_l"]
    D = inputs["D"]
    L = inputs["L"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    n_perp = inputs["n_perp"]
    Gamma_par = inputs["Gamma_par"]
    Gamma_perp = inputs["Gamma_perp"]
    eps_1 = inputs["eps_1"]
    eps_par = inputs["eps_par"]
    eps_perp = inputs["eps_perp"]
    v_l = inputs["v_l"]
    n_quad = int(settings.get("n_quad", 12))

    D_eff = D + kappa_inv * (_EULER_GAMMA + math.log(Gamma_perp))
    n_par = L / delta
    h = kappa_inv / D_eff

    # Gauss-Legendre nodes for spherical-cap integration
    xu, wu = leggauss(n_quad)
    xp, wp = leggauss(n_quad * 2)
    tu = 0.5 * (xu + 1.0)
    wtu = 0.5 * wu
    phis = 0.5 * math.pi * (xp + 1.0)
    wtp = 0.5 * wp
    W3 = wtu[:, None, None] * wtu[None, :, None] * wtp[None, None, :]

    gs_iso = np.linspace(D / L, 0.5 * math.pi, 2000)
    sg_iso = np.sin(gs_iso)

    def B_corr_array(g, u_par_tot, u_perp):
        small = g <= D / L
        out = np.empty_like(g)
        if np.any(small):
            Y = Gamma_par * (L / D) - abs(u_par_tot)
            br = 0.5 * Y if abs(Y) < 1e-8 else 1.0 + (math.exp(-Y) - 1.0) / Y
            out[small] = (
                8.0 * (D ** 2 + (kappa_inv ** 2) * (math.log(Gamma_par * L / D) ** 2)) / (D_eff * L)
                + 16.0 * (D * delta / (D_eff * L)) * br
            )
        if np.any(~small):
            sg = np.sin(g[~small])
            out[~small] = (
                (8.0 / math.pi) * sg * (1.0 - h * np.log(sg))
                + (8.0 * delta / (math.pi * D_eff)) * sg * (1.0 - np.exp(-Gamma_perp / sg + abs(u_perp)))
            )
        return out

    def eval_u(mu_l):
        u_par_tot = -n_par * _log1pexp(mu_l - 2.0 * eps_par) + 2.0 * n_par * _log1pexp(mu_l - eps_1)
        u_perp = -n_perp * _log1pexp(mu_l - 2.0 * eps_perp) + 2.0 * n_perp * _log1pexp(mu_l - eps_1)
        return u_par_tot, u_perp

    # Estimate mu_l corresponding to c_l_target in the low-c_l regime (I-N coexistence):
    # At low c_l (Fig. 3a, eps_perp=-7), free + 1-rod linkers dominate over 2-body junctions.
    def cl_approx(mu_l, cr_ref=4.0):
        return (L ** 2 * D_eff / v_l) * math.exp(mu_l) + n_par * cr_ref * _sigmoid(mu_l - eps_1)

    if c_l_target <= 0.0:
        mu_l = -40.0
    else:
        mu_l = brentq(lambda m: cl_approx(m) - c_l_target, -40.0, eps_1 + 5.0, xtol=1e-10)

    u_par_tot, u_perp = eval_u(mu_l)
    # Include cap contribution for w_I
    Y0 = Gamma_par * (L / D) - abs(u_par_tot)
    br0 = 0.5 * Y0 if abs(Y0) < 1e-8 else 1.0 + (math.exp(-Y0) - 1.0) / Y0
    B0 = (
        8.0 * (D ** 2 + (kappa_inv ** 2) * (math.log(Gamma_par * L / D) ** 2)) / (D_eff * L)
        + 16.0 * (D * delta / (D_eff * L)) * br0
    )
    w_I = (1.0 - math.cos(D / L)) * B0 + float(
        np.trapezoid(B_corr_array(gs_iso, u_par_tot, u_perp) * sg_iso, gs_iso)
    )

    def w_N(alpha):
        ca = math.cos(alpha)
        u1 = ca + (1.0 - ca) * tu
        s1 = np.sqrt(np.maximum(0.0, 1.0 - u1 ** 2))
        cg = np.abs(
            u1[:, None, None] * u1[None, :, None]
            + s1[:, None, None] * s1[None, :, None] * np.cos(phis)[None, None, :]
        )
        g = np.arccos(np.clip(cg, 0.0, 1.0))
        return float(np.sum(B_corr_array(g, u_par_tot, u_perp) * W3))

    def opt_N(cr):
        def obj(la):
            a = math.exp(la)
            sig = math.log(1.0 / (1.0 - math.cos(a)))
            return cr * sig + 0.5 * (cr ** 2) * w_N(a)

        res = minimize_scalar(obj, bounds=(math.log(1e-3), math.log(0.5 * math.pi - 1e-3)), method="bounded")
        a = math.exp(res.x)
        sig = math.log(1.0 / (1.0 - math.cos(a)))
        w = w_N(a)
        mu = math.log(cr) + sig + cr * w
        P = cr + 0.5 * (cr ** 2) * w
        return mu, P

    def diff_IN(cN):
        muN, PN = opt_N(cN)
        cI = (-1.0 + math.sqrt(1.0 + 2.0 * w_I * PN)) / w_I
        muI = math.log(cI) + cI * w_I
        return muN - muI

    c_N = brentq(diff_IN, 2.0, 30.0, xtol=1e-7)
    _, PN = opt_N(c_N)
    c_I = (-1.0 + math.sqrt(1.0 + 2.0 * w_I * PN)) / w_I
    return {"c_I": float(c_I), "c_N": float(c_N)}
