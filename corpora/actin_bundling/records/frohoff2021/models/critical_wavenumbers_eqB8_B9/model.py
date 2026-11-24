import math


def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    ell = inputs["ell"]

    f1_pp = a + 3.0 * p1 * p1
    f2_pp = a + a_Delta + 3.0 * p2 * p2
    Delta = alpha * alpha - rho * rho

    inner = (f2_pp - kappa * f1_pp) ** 2 - 4.0 * kappa * Delta
    q_plus = q_minus = k_plus = k_minus = None
    if inner >= 0.0 and kappa > 0.0:
        s = math.sqrt(inner)
        qp2 = (-f2_pp - kappa * f1_pp + s) / (2.0 * kappa)
        qm2 = (-f2_pp - kappa * f1_pp - s) / (2.0 * kappa)
        if qp2 >= 0.0:
            q_plus = math.sqrt(qp2)
            k_plus = ell * q_plus
        if qm2 >= 0.0:
            q_minus = math.sqrt(qm2)
            k_minus = ell * q_minus

    qo2 = -(2.0 * a + a_Delta + 3.0 * (p1 * p1 + p2 * p2)) / (1.0 + kappa)
    q_o = math.sqrt(qo2) if qo2 >= 0.0 else None
    k_o = ell * q_o if q_o is not None else None

    return {
        "q_plus": q_plus,
        "q_minus": q_minus,
        "k_plus": k_plus,
        "k_minus": k_minus,
        "q_o": q_o,
        "k_o": k_o,
    }


def symbolic():
    import sympy as sp

    a, a_Delta, kappa, phi1_bar, phi2_bar, rho, alpha, ell = sp.symbols(
        "a a_Delta kappa phi1_bar phi2_bar rho alpha ell", real=True
    )
    f1_pp = a + 3 * phi1_bar**2
    f2_pp = a + a_Delta + 3 * phi2_bar**2
    Delta = alpha**2 - rho**2
    s = sp.sqrt((f2_pp - kappa * f1_pp) ** 2 - 4 * kappa * Delta)
    q_plus = sp.sqrt((-f2_pp - kappa * f1_pp + s) / (2 * kappa))
    q_minus = sp.sqrt((-f2_pp - kappa * f1_pp - s) / (2 * kappa))
    q_o = sp.sqrt(-(2 * a + a_Delta + 3 * (phi1_bar**2 + phi2_bar**2)) / (1 + kappa))
    return {
        "q_plus": q_plus,
        "q_minus": q_minus,
        "k_plus": ell * q_plus,
        "k_minus": ell * q_minus,
        "q_o": q_o,
        "k_o": ell * q_o,
    }
