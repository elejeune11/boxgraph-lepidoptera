import math


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta"])
    arg = (k * k - d * d - aA * aB) / aB
    if arg < 0.0:
        return {
            "rho_A_s": None,
            "rho_B_s": None,
            "theta_s": math.pi,
            "amp_ratio_s": aB / (k + d),
            "v_s": 0.0,
        }
    rho_A = math.sqrt(arg)
    rho_B = (k + d) * rho_A / aB
    return {
        "rho_A_s": rho_A,
        "rho_B_s": rho_B,
        "theta_s": math.pi,
        "amp_ratio_s": aB / (k + d),
        "v_s": 0.0,
    }


def symbolic():
    import sympy as sp

    alpha_A, alpha_B, kappa, delta = sp.symbols("alpha_A alpha_B kappa delta", real=True)
    rho_A = sp.sqrt((kappa**2 - delta**2 - alpha_A * alpha_B) / alpha_B)
    rho_B = (kappa + delta) * rho_A / alpha_B
    return {
        "rho_A_s": rho_A,
        "rho_B_s": rho_B,
        "theta_s": sp.pi,
        "amp_ratio_s": alpha_B / (kappa + delta),
        "v_s": sp.Integer(0),
    }
