import math

def run(inputs, settings):
    kd = inputs["k_d"]
    skd = math.sqrt(kd)
    b = inputs["b"]
    beta = inputs["beta"]
    lo = inputs["lambda_o"]
    kap = inputs["kappa"]
    q0 = inputs["q_0"]

    delta = kap * beta / (16.0 * b * q0 * q0)
    tau = kap * (lo / (8.0 * b * q0) + q0)
    op = 1.0 + tau
    nu0_sq = skd / 2.0
    lam0 = (1.0 + 2.0 * skd) ** 2
    nu_sq = nu0_sq * (1.0 + 0.125 * (((1.0 + 2.0 * skd) ** 2) / op) * delta)
    lam = (lam0 / op) * (1.0 - 0.5 * skd * ((1.0 + 2.0 * skd) / op) * delta)
    nu = math.sqrt(nu_sq)
    return {
        "delta": delta,
        "tau": tau,
        "nu_crit_sq": nu_sq,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2.0 * math.pi / nu,
    }

def symbolic():
    import sympy as sp
    k_d, b, beta, lambda_o, kappa, q_0 = sp.symbols("k_d b beta lambda_o kappa q_0", real=True)
    skd = sp.sqrt(k_d)
    delta = kappa * beta / (16 * b * q_0**2)
    tau = kappa * (lambda_o / (8 * b * q_0) + q_0)
    op = 1 + tau
    nu0_sq = skd / 2
    lam0 = (1 + 2 * skd) ** 2
    nu_sq = nu0_sq * (1 + sp.Rational(1, 8) * ((1 + 2 * skd) ** 2 / op) * delta)
    lam = (lam0 / op) * (1 - skd / 2 * ((1 + 2 * skd) / op) * delta)
    nu = sp.sqrt(nu_sq)
    return {
        "delta": delta,
        "tau": tau,
        "nu_crit_sq": nu_sq,
        "nu_crit": nu,
        "lambda_crit": lam,
        "ell_pattern": 2 * sp.pi / nu,
    }
