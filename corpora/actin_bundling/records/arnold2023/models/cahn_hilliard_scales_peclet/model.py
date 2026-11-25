import math


def run(inputs, settings):
    M = inputs["M"]
    kappa = inputs["kappa"]
    gamma_dot = inputs["gamma_dot"]
    return {
        "l_c": math.sqrt(kappa),
        "t_c": kappa / M,
        "Pe": gamma_dot * kappa / M,
    }


def symbolic():
    import sympy as sp
    M, kappa, gamma_dot = sp.symbols("M kappa gamma_dot", positive=True)
    return {
        "l_c": sp.sqrt(kappa),
        "t_c": kappa / M,
        "Pe": gamma_dot * kappa / M,
    }
