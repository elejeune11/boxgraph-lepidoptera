def run(inputs, settings):
    return {"y": inputs["R_tilde"] / inputs["R"]}


def symbolic():
    import sympy as sp

    R_tilde, R = sp.symbols("R_tilde R", positive=True)
    return {"y": R_tilde / R}
