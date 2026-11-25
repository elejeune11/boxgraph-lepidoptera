def run(inputs, settings):
    return {"nu_bar": (inputs["eta_s"] + inputs["eta_b"]) / (4.0 * inputs["eta_s"])}


def symbolic():
    import sympy as sp

    eta_s, eta_b = sp.symbols("eta_s eta_b", positive=True)
    return {"nu_bar": (eta_s + eta_b) / (4 * eta_s)}
