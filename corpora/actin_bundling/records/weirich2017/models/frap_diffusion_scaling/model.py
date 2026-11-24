def run(inputs, settings):
    return {"D": inputs["L_b"] ** 2 / inputs["tau_R"]}


def symbolic():
    import sympy as sp

    L_b, tau_R = sp.symbols("L_b tau_R", positive=True)
    return {"D": L_b**2 / tau_R}
