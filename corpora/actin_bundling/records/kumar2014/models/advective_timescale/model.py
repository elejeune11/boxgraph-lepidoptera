def run(inputs, settings):
    return {"tau_a": inputs["ell"] / inputs["U"]}


def symbolic():
    import sympy as sp

    ell, U = sp.symbols("ell U", positive=True)
    return {"tau_a": ell / U}
