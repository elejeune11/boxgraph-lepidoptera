import math


def run(inputs, settings):
    a = inputs["a"]
    kappa = inputs["kappa"]
    q_n = inputs["q_n"]

    phi2_q_sq = 0.5 * kappa * q_n * q_n
    phi2_q = math.sqrt(phi2_q_sq) if phi2_q_sq >= 0.0 else None
    B2 = 9.0 * phi2_q_sq if phi2_q_sq >= 0.0 else None
    phi2_a = math.sqrt(-a / 5.0) if a <= 0.0 else None
    return {
        "phi2_bar_crit_q": phi2_q,
        "B2_crit": B2,
        "phi2_bar_crit_a": phi2_a,
    }


def symbolic():
    import sympy as sp

    a, kappa, q_n = sp.symbols("a kappa q_n", real=True)
    phi2_q_sq = sp.Rational(1, 2) * kappa * q_n**2
    return {
        "phi2_bar_crit_q": sp.sqrt(phi2_q_sq),
        "B2_crit": 9 * phi2_q_sq,
        "phi2_bar_crit_a": sp.sqrt(-a / 5),
    }
