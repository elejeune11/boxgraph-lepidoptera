import math

def run(inputs, settings):
    gamma = inputs["gamma"]
    D_si = inputs["D"] * 1e-12
    kd = inputs["k_d"]
    eta = inputs["eta"]
    lam0 = (math.sqrt(gamma * D_si) + 2.0 * math.sqrt(kd * eta)) ** 2
    nu0_sq = math.sqrt(kd * gamma / (4.0 * eta * D_si)) * 1e-12
    nu0 = math.sqrt(nu0_sq)
    return {
        "lambda_crit_0": lam0,
        "nu_crit_0_sq": nu0_sq,
        "nu_crit_0": nu0,
        "ell_pattern_0": 2.0 * math.pi / nu0,
    }

def symbolic():
    import sympy as sp
    gamma, D, k_d, eta = sp.symbols("gamma D k_d eta", positive=True)
    D_si = D * sp.Rational(1, 10**12)
    lam0 = (sp.sqrt(gamma * D_si) + 2 * sp.sqrt(k_d * eta)) ** 2
    nu0_sq = sp.sqrt(k_d * gamma / (4 * eta * D_si)) * sp.Rational(1, 10**12)
    nu0 = sp.sqrt(nu0_sq)
    return {
        "lambda_crit_0": lam0,
        "nu_crit_0_sq": nu0_sq,
        "nu_crit_0": nu0,
        "ell_pattern_0": 2 * sp.pi / nu0,
    }
