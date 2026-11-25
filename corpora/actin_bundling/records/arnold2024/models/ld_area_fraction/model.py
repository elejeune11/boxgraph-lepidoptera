def run(inputs, settings):
    phi_Lo = inputs["phi_Lo"]
    phi_Ld = 100.0 - phi_Lo
    delta_phi_total = phi_Lo
    delta_phi_rel = 100.0 * phi_Lo / phi_Ld
    return {
        "phi_Ld": phi_Ld,
        "delta_phi_total": delta_phi_total,
        "delta_phi_rel": delta_phi_rel,
    }


def symbolic():
    import sympy as sp

    phi_Lo = sp.symbols("phi_Lo", positive=True)
    phi_Ld = sp.Integer(100) - phi_Lo
    return {
        "phi_Ld": phi_Ld,
        "delta_phi_total": phi_Lo,
        "delta_phi_rel": sp.Integer(100) * phi_Lo / phi_Ld,
    }
