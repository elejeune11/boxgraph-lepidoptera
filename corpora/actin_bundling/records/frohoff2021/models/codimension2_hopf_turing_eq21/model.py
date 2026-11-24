import math


def run(inputs, settings):
    Q = inputs["Q"]
    kappa = inputs["kappa"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    Delta = alpha * alpha - rho * rho
    if Delta < 0.0 or kappa < 0.0:
        return {"f2_pp_cd2": None, "f1_pp_cd2": None}
    sign = 1.0 if kappa >= 1.0 / Q else -1.0
    f2_cd2 = 2.0 * sign * math.sqrt(kappa * Delta) / (1.0 + Q * kappa)
    f1_cd2 = -Q * f2_cd2
    return {"f2_pp_cd2": f2_cd2, "f1_pp_cd2": f1_cd2}


def symbolic():
    import sympy as sp

    Q, kappa, rho, alpha = sp.symbols("Q kappa rho alpha", real=True)
    Delta = alpha**2 - rho**2
    sign = sp.Piecewise((1, kappa * Q >= 1), (-1, True))
    f2_cd2 = 2 * sign * sp.sqrt(kappa * Delta) / (1 + Q * kappa)
    return {
        "f2_pp_cd2": f2_cd2,
        "f1_pp_cd2": -Q * f2_cd2,
    }
