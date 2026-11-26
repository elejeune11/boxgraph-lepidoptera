import math

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
    L = inputs["L"]
    D = inputs["D"]
    delta = inputs["delta"]
    kappa_inv = inputs["kappa_inv"]
    Gamma_perp = inputs["Gamma_perp"]
    v_l = inputs["v_l"]
    c_r = inputs["c_r"]
    mu_l = inputs["mu_l"]
    eps_1 = inputs["eps_1"]

    D_eff = D + kappa_inv * (_EULER_GAMMA + math.log(Gamma_perp))
    n_par = L / delta
    free_term = (L ** 2 * D_eff / v_l) * math.exp(mu_l)
    f_linkers = -free_term - n_par * c_r * _log1pexp(mu_l - eps_1)
    c_l_non_junction = free_term + n_par * c_r * _sigmoid(mu_l - eps_1)
    return {
        "f_linkers": float(f_linkers),
        "c_l_non_junction": float(c_l_non_junction),
    }


def symbolic():
    import sympy as sp

    L, D, delta, kappa_inv, Gamma_perp, v_l, c_r = sp.symbols(
        "L D delta kappa_inv Gamma_perp v_l c_r", positive=True
    )
    mu_l, eps_1 = sp.symbols("mu_l eps_1", real=True)
    D_eff = D + kappa_inv * (sp.EulerGamma + sp.log(Gamma_perp))
    n_par = L / delta
    free_term = L ** 2 * D_eff / v_l * sp.exp(mu_l)
    f_linkers = -free_term - n_par * c_r * sp.log(1 + sp.exp(mu_l - eps_1))
    c_l_non_junction = free_term + n_par * c_r * sp.exp(mu_l - eps_1) / (1 + sp.exp(mu_l - eps_1))
    return {"f_linkers": f_linkers, "c_l_non_junction": c_l_non_junction}
