import math


def run(inputs, settings):
    rho_star = float(inputs["rho_star"])
    GR = float(inputs["Gamma_R"])
    GA = float(inputs["Gamma_A"])
    lA = float(inputs["lambda_A"])
    # From beta V(0) = 2*pi*GR*lR^2 - pi*GA*lA^2 < 0 and rho* * (-beta V(0)) * b(0) = 1 (with b(0) = 1)
    num = math.pi * GA * (lA * lA) - 1.0 / rho_star
    if num <= 0.0:
        return {"lambda_R": None}
    return {"lambda_R": math.sqrt(num / (2.0 * math.pi * GR))}


def symbolic():
    import sympy as sp
    rho_star, Gamma_R, Gamma_A, lambda_A = sp.symbols(
        "rho_star Gamma_R Gamma_A lambda_A", positive=True
    )
    return {"lambda_R": sp.sqrt((sp.pi * Gamma_A * lambda_A**2 - 1 / rho_star) / (2 * sp.pi * Gamma_R))}
