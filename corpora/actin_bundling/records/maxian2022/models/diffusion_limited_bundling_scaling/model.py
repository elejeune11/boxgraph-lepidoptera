def run(inputs, settings):
    # Unknown prefactor set to 1 (in s / um^2); BoxGraph fits the prefactor against curve data.
    return {"tau_c": inputs["ell_m"] ** 2}


def symbolic():
    import sympy as sp

    ell_m = sp.Symbol("ell_m", positive=True)
    return {"tau_c": ell_m**2}
