import math


def run(inputs, settings):
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    Delta = alpha * alpha - rho * rho
    a_H = -0.5 * (a_Delta + 3.0 * (p1 * p1 + p2 * p2))
    diff0 = 3.0 * (p1 * p1 - p2 * p2) - a_Delta
    disc0 = diff0 * diff0 - 4.0 * Delta
    a_CH = a_H + 0.5 * math.sqrt(disc0) if disc0 >= 0.0 else None

    if Delta > 0.0 and kappa > 0.0 and abs(kappa - 1.0) > 1e-15:
        sign = -1.0 if kappa > 1.0 else 1.0
        a_T = (a_Delta + 3.0 * (p2 * p2 - kappa * p1 * p1) + sign * 2.0 * math.sqrt(kappa * Delta)) / (kappa - 1.0)
    else:
        a_T = None

    return {"a_H": a_H, "a_CH": a_CH, "a_T": a_T}


def symbolic():
    import sympy as sp

    a_Delta, kappa, phi1_bar, phi2_bar, rho, alpha = sp.symbols(
        "a_Delta kappa phi1_bar phi2_bar rho alpha", real=True
    )
    Delta = alpha**2 - rho**2
    a_H = -sp.Rational(1, 2) * (a_Delta + 3 * (phi1_bar**2 + phi2_bar**2))
    diff0 = 3 * (phi1_bar**2 - phi2_bar**2) - a_Delta
    a_CH = a_H + sp.Rational(1, 2) * sp.sqrt(diff0**2 - 4 * Delta)
    sign = sp.Piecewise((-1, kappa > 1), (1, True))
    a_T = (a_Delta + 3 * (phi2_bar**2 - kappa * phi1_bar**2) + sign * 2 * sp.sqrt(kappa * Delta)) / (kappa - 1)
    return {"a_H": a_H, "a_CH": a_CH, "a_T": a_T}
