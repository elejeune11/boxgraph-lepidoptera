import math


def run(inputs, settings):
    f2_pp = inputs["f2_pp"]
    Q = inputs["Q"]
    kappa = inputs["kappa"]
    Delta = inputs["Delta"]

    # Eq. (14): Hopf threshold f_1''^H = -Q * f_2'' for Delta > Q * f_2''^2
    f1_pp_H = -Q * f2_pp if Delta > Q * f2_pp * f2_pp else None

    # Eq. (19): Turing threshold (relevant stability boundary branch: - for kappa > 1/Q, + for kappa < 1/Q)
    f1_pp_T = None
    if Delta > 0.0 and kappa > 0.0 and abs(Q * kappa - 1.0) > 1e-15:
        sign = 1.0 if kappa > 1.0 / Q else -1.0
        skd = math.sqrt(kappa * Delta)
        f1_pp_T = (f2_pp - 2.0 * sign * skd) / kappa

    # Eq. (20): CH threshold f_1''^CH = -Delta / f_2''
    f1_pp_CH = -Delta / f2_pp if abs(f2_pp) > 1e-15 else None

    return {
        "f1_pp_H": f1_pp_H,
        "f1_pp_T": f1_pp_T,
        "f1_pp_CH": f1_pp_CH,
    }


def symbolic():
    import sympy as sp

    f2_pp, Q, kappa, Delta = sp.symbols("f2_pp Q kappa Delta", real=True)
    sign = sp.Piecewise((1, kappa * Q > 1), (-1, True))
    return {
        "f1_pp_H": -Q * f2_pp,
        "f1_pp_T": (f2_pp - 2 * sign * sp.sqrt(kappa * Delta)) / kappa,
        "f1_pp_CH": -Delta / f2_pp,
    }
