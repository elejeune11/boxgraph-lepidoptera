import math

def run(inputs, settings):
    kd = inputs["k_d"]
    delta = inputs["kappa"] * inputs["beta"] / (2.0 * inputs["L"])
    x = (-2.0 * kd * delta + math.sqrt(kd * (1.0 + delta * (1.0 - 4.0 * kd)))) / (2.0 * (1.0 + delta))
    lam = (1.0 + 4.0 * x) * (kd + x) / ((1.0 + delta) * x + kd * delta)
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
    k_d, L, beta, kappa = sp.symbols("k_d L beta kappa", real=True)
    delta = kappa * beta / (2 * L)
    x = (-2 * k_d * delta + sp.sqrt(k_d * (1 + delta * (1 - 4 * k_d)))) / (2 * (1 + delta))
    lam = (1 + 4 * x) * (k_d + x) / ((1 + delta) * x + k_d * delta)
    nu = sp.sqrt(x)
    return {
        "delta": delta,
        "nu_crit_sq": x,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2 * sp.pi / nu,
    }
