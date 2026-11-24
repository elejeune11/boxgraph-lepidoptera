def run(inputs, settings):
    return {"alpha_visc": inputs["eta_drop"] / inputs["eta_water"]}


def symbolic():
    import sympy as sp

    eta_drop, eta_water = sp.symbols("eta_drop eta_water", positive=True)
    return {"alpha_visc": eta_drop / eta_water}
