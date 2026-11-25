import math

def run(inputs, settings):
    kd = inputs["k_d"]
    skd = math.sqrt(kd)
    delta = inputs["kappa"] * inputs["beta"] / (2.0 * inputs["c_0"])
    nu0_sq = skd / 2.0
    lam0 = (1.0 + 2.0 * skd) ** 2
    nu_sq = nu0_sq * (1.0 + 0.125 * ((1.0 + 2.0 * skd) ** 2) * delta)
    lam = lam0 * (1.0 - 0.5 * skd * (1.0 + 2.0 * skd) * delta)
    nu = math.sqrt(nu_sq)
    return {
        "delta": delta,
        "nu_crit_sq": nu_sq,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2.0 * math.pi / nu,
    }

def symbolic():
    import sympy as sp
    k_d, c_0, beta, kappa = sp.symbols("k_d c_0 beta kappa", real=True)
    skd = sp.sqrt(k_d)
    delta = kappa * beta / (2 * c_0)
    nu0_sq = skd / 2
    lam0 = (1 + 2 * skd) ** 2
    nu_sq = nu0_sq * (1 + sp.Rational(1, 8) * (1 + 2 * skd) ** 2 * delta)
    lam = lam0 * (1 - skd / 2 * (1 + 2 * skd) * delta)
    nu = sp.sqrt(nu_sq)
    return {
        "delta": delta,
        "nu_crit_sq": nu_sq,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2 * sp.pi / nu,
    }
