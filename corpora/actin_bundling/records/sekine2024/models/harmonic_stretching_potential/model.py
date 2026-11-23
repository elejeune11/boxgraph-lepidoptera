"""Harmonic stretching potential U^harm(r; kappa, L_0) = (kappa / (2 L_0)) (r - L_0)^2 (Eq. S5)."""


def run(inputs, settings):
    r = inputs["r"]
    kappa = inputs["kappa"]
    L0 = inputs["L0"]
    k_eff = kappa / L0
    U_harm = 0.5 * k_eff * (r - L0) ** 2
    f_str = -k_eff * (r - L0)
    return {"U_harm": U_harm, "f_str": f_str, "k_eff": k_eff}


def symbolic():
    import sympy as sp

    r, kappa, L0 = sp.symbols("r kappa L0", positive=True)
    return {
        "U_harm": (kappa / (2 * L0)) * (r - L0) ** 2,
        "f_str": -(kappa / L0) * (r - L0),
        "k_eff": kappa / L0,
    }
