"""Harmonic stretching potential for motor/crosslinker half-arm rest length L_0/2 (Eqs. S6, S8)."""


def run(inputs, settings):
    r = inputs["r"]
    kappa = inputs["kappa"]
    half_L0 = 0.5 * inputs["L0"]
    k_eff = kappa / half_L0
    U_harm = 0.5 * k_eff * (r - half_L0) ** 2
    f_str = -k_eff * (r - half_L0)
    return {"U_harm": U_harm, "f_str": f_str, "k_eff": k_eff}


def symbolic():
    import sympy as sp

    r, kappa, L0 = sp.symbols("r kappa L0", positive=True)
    return {
        "U_harm": (kappa / L0) * (r - L0 / 2) ** 2,
        "f_str": -(2 * kappa / L0) * (r - L0 / 2),
        "k_eff": 2 * kappa / L0,
    }
