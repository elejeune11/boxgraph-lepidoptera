def run(inputs, settings):
    # (Pa*s) * (um/s) = (N/m^2 * s) * (1e-6 m/s) = 1e-6 N/m = 1e3 nN/m
    return {"gamma_I": inputs["eta"] * inputs["v_cap"] * 1e3}


def symbolic():
    import sympy as sp

    eta, v_cap = sp.symbols("eta v_cap", positive=True)
    return {"gamma_I": sp.Integer(1000) * eta * v_cap}
