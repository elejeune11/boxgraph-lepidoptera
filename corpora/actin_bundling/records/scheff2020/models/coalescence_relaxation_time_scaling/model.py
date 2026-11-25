def run(inputs, settings):
    return {"tau": inputs["eta_over_gamma"] * inputs["delta_l"]}


def symbolic():
    import sympy as sp

    eta_over_gamma, delta_l = sp.symbols("eta_over_gamma delta_l", positive=True)
    return {"tau": eta_over_gamma * delta_l}
