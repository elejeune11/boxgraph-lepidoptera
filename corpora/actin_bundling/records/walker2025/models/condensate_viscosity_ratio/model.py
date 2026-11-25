def run(inputs, settings):
    # 1 pN*s/um^2 = 1e-12 N * s / 1e-12 m^2 = 1 Pa*s exactly
    return {"ratio": float(inputs["eta_cond"]) / float(inputs["eta_water"])}


def symbolic():
    import sympy as sp

    eta_cond, eta_water = sp.symbols("eta_cond eta_water", positive=True)
    return {"ratio": eta_cond / eta_water}
