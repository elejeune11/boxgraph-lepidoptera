def run(inputs, settings):
    return {"gamma_beta": inputs["w"] * inputs["psi"] / 6.0}


def symbolic():
    import sympy as sp

    w, psi = sp.symbols("w psi", positive=True)
    return {"gamma_beta": w * psi / 6}
