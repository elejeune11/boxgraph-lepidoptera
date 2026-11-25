import math

def run(inputs, settings):
    gamma_D = inputs["gamma"] * inputs["D"] * 1e-12
    lam0 = gamma_D * (1.0 + 2.0 * math.sqrt(inputs["ell_s"] / inputs["ell_D"])) ** 2
    nu0_sq = 1.0 / (2.0 * inputs["ell_s"] * inputs["ell_D"])
    return {"lambda_crit_0": lam0, "nu_crit_0_sq": nu0_sq}

def symbolic():
    import sympy as sp
    gamma, D, ell_s, ell_D = sp.symbols("gamma D ell_s ell_D", positive=True)
    gamma_D = gamma * D * sp.Rational(1, 10**12)
    return {
        "lambda_crit_0": gamma_D * (1 + 2 * sp.sqrt(ell_s / ell_D)) ** 2,
        "nu_crit_0_sq": 1 / (2 * ell_s * ell_D),
    }
