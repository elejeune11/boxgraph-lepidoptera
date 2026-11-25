import math

def run(inputs, settings):
    gamma = inputs["gamma"]
    D_si = inputs["D"] * 1e-12
    kd = inputs["k_d"]
    eta = inputs["eta"]
    kap = inputs["kappa"]
    beta = inputs["beta"]
    c0 = inputs["c_0"]

    ell_s = math.sqrt(eta / gamma) * 1e6
    ell_D = math.sqrt(inputs["D"] / kd)
    r = ell_s / ell_D
    delta = gamma * D_si * kap * beta / (2.0 * eta * c0)
    lam0 = (math.sqrt(gamma * D_si) + 2.0 * math.sqrt(kd * eta)) ** 2
    lam_c = lam0 * (1.0 - 0.5 * r * (1.0 + 2.0 * r) * delta)
    nu0_sq = 1.0 / (2.0 * ell_s * ell_D)
    nu_sq = nu0_sq * (1.0 + 0.125 * ((1.0 + 2.0 * r) ** 2) * delta)
    nu = math.sqrt(nu_sq)
    return {
        "delta": delta,
        "lambda_crit": lam_c,
        "nu_crit_sq": nu_sq,
        "nu_crit": nu,
        "ell_pattern": 2.0 * math.pi / nu,
    }

def symbolic():
    import sympy as sp
    gamma, D, k_d, eta, kappa, beta, c_0 = sp.symbols("gamma D k_d eta kappa beta c_0", positive=True)
    D_si = D * sp.Rational(1, 10**12)
    ell_s = sp.sqrt(eta / gamma) * 10**6
    ell_D = sp.sqrt(D / k_d)
    r = ell_s / ell_D
    delta = gamma * D_si * kappa * beta / (2 * eta * c_0)
    lam0 = (sp.sqrt(gamma * D_si) + 2 * sp.sqrt(k_d * eta)) ** 2
    lam_c = lam0 * (1 - sp.Rational(1, 2) * r * (1 + 2 * r) * delta)
    nu0_sq = 1 / (2 * ell_s * ell_D)
    nu_sq = nu0_sq * (1 + sp.Rational(1, 8) * (1 + 2 * r) ** 2 * delta)
    nu = sp.sqrt(nu_sq)
    return {
        "delta": delta,
        "lambda_crit": lam_c,
        "nu_crit_sq": nu_sq,
        "nu_crit": nu,
        "ell_pattern": 2 * sp.pi / nu,
    }
