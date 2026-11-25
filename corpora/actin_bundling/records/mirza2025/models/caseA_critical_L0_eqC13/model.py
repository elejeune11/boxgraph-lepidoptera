import math

def run(inputs, settings):
    kd = inputs["k_d"]
    delta = inputs["kappa"] * inputs["beta"] / (2.0 * inputs["c_0"])
    x = (kd * delta + math.sqrt(kd * (4.0 + delta * (4.0 * kd - 1.0)))) / (4.0 - delta)
    lam = (1.0 + 4.0 * x) * (kd + x) / (delta * x * x + (1.0 + kd * delta) * x)
    nu = math.sqrt(x)
    return {
        "delta": delta,
        "nu_crit_sq": x,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2.0 * math.pi / nu,
    }

def symbolic():
    import sympy as sp
    k_d, c_0, beta, kappa = sp.symbols("k_d c_0 beta kappa", real=True)
    delta = kappa * beta / (2 * c_0)
    x = (k_d * delta + sp.sqrt(k_d * (4 + delta * (4 * k_d - 1)))) / (4 - delta)
    lam = (1 + 4 * x) * (k_d + x) / (delta * x**2 + (1 + k_d * delta) * x)
    nu = sp.sqrt(x)
    return {
        "delta": delta,
        "nu_crit_sq": x,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2 * sp.pi / nu,
    }
