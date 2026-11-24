import math


def run(inputs, settings):
    f2_pp = inputs["f2_pp"]
    Q = inputs["Q"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    Delta = alpha * alpha - rho * rho
    f1_pp_H = -Q * f2_pp
    rad = -Q * f2_pp * f2_pp + Delta
    omega_tilde_c = math.sqrt(Q * rad) if rad >= 0.0 else None
    f2_pp_H_end = math.sqrt(Delta / Q) if Delta >= 0.0 else None
    return {
        "f1_pp_H": f1_pp_H,
        "omega_tilde_c": omega_tilde_c,
        "f2_pp_H_end": f2_pp_H_end,
    }


def symbolic():
    import sympy as sp

    f2_pp, Q, rho, alpha = sp.symbols("f2_pp Q rho alpha", real=True)
    Delta = alpha**2 - rho**2
    return {
        "f1_pp_H": -Q * f2_pp,
        "omega_tilde_c": sp.sqrt(Q) * sp.sqrt(-Q * f2_pp**2 + Delta),
        "f2_pp_H_end": sp.sqrt(Delta / Q),
    }
