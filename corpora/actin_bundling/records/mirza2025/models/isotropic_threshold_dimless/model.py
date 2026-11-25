import math

def run(inputs, settings):
    kd = inputs["k_d"]
    nu_sq = math.sqrt(kd) / 2.0
    nu = math.sqrt(nu_sq)
    lam = (1.0 + 2.0 * math.sqrt(kd)) ** 2
    return {
        "nu_crit_0_sq": nu_sq,
        "nu_crit_0": nu,
        "lambda_crit_0": lam,
        "ell_pattern_0": 2.0 * math.pi / nu,
    }

def symbolic():
    import sympy as sp
    k_d = sp.symbols("k_d", positive=True)
    nu_sq = sp.sqrt(k_d) / 2
    nu = sp.sqrt(nu_sq)
    lam = (1 + 2 * sp.sqrt(k_d)) ** 2
    return {
        "nu_crit_0_sq": nu_sq,
        "nu_crit_0": nu,
        "lambda_crit_0": lam,
        "ell_pattern_0": 2 * sp.pi / nu,
    }
