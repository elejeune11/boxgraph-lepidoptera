import math

def run(inputs, settings):
    kd = inputs["k_d"]
    b = inputs["b"]
    beta = inputs["beta"]
    lo = inputs["lambda_o"]
    kap = inputs["kappa"]
    q0 = inputs["q_0"]

    delta = kap * beta / (16.0 * b * q0 * q0)
    tau = kap * (lo / (8.0 * b * q0) + q0)
    denom = 4.0 + 4.0 * tau - delta
    x = (delta * kd + math.sqrt((delta * kd) ** 2 + kd * denom * (1.0 + tau + delta * kd))) / denom
    lam = (1.0 + 4.0 * x) * (x + kd) / (delta * x * x + (1.0 + tau + delta * kd) * x)
    nu = math.sqrt(x)
    return {
        "delta": delta,
        "tau": tau,
        "nu_crit_sq": x,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2.0 * math.pi / nu,
    }

def symbolic():
    import sympy as sp
    k_d, b, beta, lambda_o, kappa, q_0 = sp.symbols("k_d b beta lambda_o kappa q_0", real=True)
    delta = kappa * beta / (16 * b * q_0**2)
    tau = kappa * (lambda_o / (8 * b * q_0) + q_0)
    denom = 4 + 4 * tau - delta
    x = (delta * k_d + sp.sqrt(delta**2 * k_d**2 + k_d * denom * (1 + tau + delta * k_d))) / denom
    lam = (1 + 4 * x) * (x + k_d) / (delta * x**2 + (1 + tau + delta * k_d) * x)
    nu = sp.sqrt(x)
    return {
        "delta": delta,
        "tau": tau,
        "nu_crit_sq": x,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2 * sp.pi / nu,
    }
