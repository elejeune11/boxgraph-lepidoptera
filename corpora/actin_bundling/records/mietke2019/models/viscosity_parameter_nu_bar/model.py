def run(inputs, settings):
    # Verbatim as printed on SI p. 15 (PDF p. 21): missing division slash.
    return {"nu_bar": (inputs["eta_s"] + inputs["eta_b"]) * (4.0 * inputs["eta_s"])}


def symbolic():
    import sympy as sp

    eta_s, eta_b = sp.symbols("eta_s eta_b", positive=True)
    return {"nu_bar": (eta_s + eta_b) * (4 * eta_s)}
