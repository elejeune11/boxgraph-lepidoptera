import math


def run(inputs, settings):
    beta = inputs["beta"]
    P = inputs["P"]
    alpha = inputs["alpha"]
    q = (inputs["omega"] * math.pi / inputs["L_over_l"]) ** 2
    Pi = P / (1.0 + inputs["nu"] * q)
    tr = -q * ((1.0 + alpha) - 2.0 * Pi)
    disc = (q**2) * (
        (1.0 - alpha) ** 2
        + 4.0 * (Pi**2)
        - 4.0 * Pi * (1.0 - alpha) * beta
    )
    if disc >= 0.0:
        raise ValueError(f"Mode k1 is not oscillatory at beta={beta}, P={P} (disc={disc})")
    im_lam = 0.5 * math.sqrt(-disc)
    return {
        "T_over_tau": 2.0 * math.pi / im_lam,
        "im_lam_k1": im_lam,
        "re_lam_k1": 0.5 * tr,
    }


def symbolic():
    import sympy as sp

    beta, P, alpha, omega, L_over_l, nu = sp.symbols(
        "beta P alpha omega L_over_l nu", positive=True
    )
    q = (omega * sp.pi / L_over_l) ** 2
    Pi = P / (1 + nu * q)
    tr = -q * ((1 + alpha) - 2 * Pi)
    neg_disc = q**2 * (4 * Pi * (1 - alpha) * beta - (1 - alpha) ** 2 - 4 * Pi**2)
    im_lam = sp.sqrt(neg_disc) / 2
    return {
        "T_over_tau": 2 * sp.pi / im_lam,
        "im_lam_k1": im_lam,
        "re_lam_k1": tr / 2,
    }
