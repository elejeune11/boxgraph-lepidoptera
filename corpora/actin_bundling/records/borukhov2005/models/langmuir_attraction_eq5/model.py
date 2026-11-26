import math


def _log1pexp(x):
    if x > 50.0:
        return x
    if x < -50.0:
        return math.exp(x)
    return math.log1p(math.exp(x))


def run(inputs, settings):
    L = inputs["L"]
    D = inputs["D"]
    delta = inputs["delta"]
    n_perp = inputs["n_perp"]
    mu_l = inputs["mu_l"]
    eps_1 = inputs["eps_1"]
    eps_par = inputs["eps_par"]
    eps_perp = inputs["eps_perp"]

    n_par = L / delta
    u_par_total = -n_par * _log1pexp(mu_l - 2.0 * eps_par) + 2.0 * n_par * _log1pexp(mu_l - eps_1)
    u_par = u_par_total * (D / L)
    u_perp = -n_perp * _log1pexp(mu_l - 2.0 * eps_perp) + 2.0 * n_perp * _log1pexp(mu_l - eps_1)
    return {
        "u_par_total": float(u_par_total),
        "u_par": float(u_par),
        "u_perp": float(u_perp),
    }


def symbolic():
    import sympy as sp

    L, D, delta, n_perp = sp.symbols("L D delta n_perp", positive=True)
    mu_l, eps_1, eps_par, eps_perp = sp.symbols("mu_l eps_1 eps_par eps_perp", real=True)
    n_par = L / delta
    u_par_total = -n_par * sp.log(1 + sp.exp(mu_l - 2 * eps_par)) + 2 * n_par * sp.log(1 + sp.exp(mu_l - eps_1))
    u_par = u_par_total * D / L
    u_perp = -n_perp * sp.log(1 + sp.exp(mu_l - 2 * eps_perp)) + 2 * n_perp * sp.log(1 + sp.exp(mu_l - eps_1))
    return {"u_par_total": u_par_total, "u_par": u_par, "u_perp": u_perp}
