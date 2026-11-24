import math


def run(inputs, settings):
    a = inputs["a"]
    f2_pp = inputs["f2_pp"]
    Q = inputs["Q"]
    kappa = inputs["kappa"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    Delta = alpha * alpha - rho * rho
    if Delta < 0.0 or kappa <= 0.0 or abs(Q * kappa - 1.0) < 1e-15:
        return {"phi1_bar_T": None, "Delta_over_kappa_max": None}

    sign = 1.0 if kappa > 1.0 / Q else -1.0
    f1_pp_T = (f2_pp - 2.0 * sign * math.sqrt(kappa * Delta)) / kappa
    rad = (f1_pp_T - a) / 3.0
    phi1_bar_T = math.sqrt(rad) if rad >= 0.0 else 0.0
    delta_over_kappa_max = 0.25 * (f2_pp / kappa - a) ** 2
    return {
        "phi1_bar_T": phi1_bar_T,
        "Delta_over_kappa_max": delta_over_kappa_max,
    }


def symbolic():
    import sympy as sp

    a, f2_pp, Q, kappa, rho, alpha = sp.symbols("a f2_pp Q kappa rho alpha", real=True)
    Delta = alpha**2 - rho**2
    sign = sp.Piecewise((1, kappa * Q > 1), (-1, True))
    f1_pp_T = (f2_pp - 2 * sign * sp.sqrt(kappa * Delta)) / kappa
    return {
        "phi1_bar_T": sp.sqrt((f1_pp_T - a) / 3),
        "Delta_over_kappa_max": sp.Rational(1, 4) * (f2_pp / kappa - a) ** 2,
    }
