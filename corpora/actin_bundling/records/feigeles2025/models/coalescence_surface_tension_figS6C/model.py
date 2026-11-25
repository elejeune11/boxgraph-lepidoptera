def run(inputs, settings):
    # (1 Pa*s) * (1 um/s) = (N*s/m^2) * (1e-6 m/s) = 1e-6 N/m = 1 uN/m
    return {"gamma": inputs["eta"] * inputs["slope"]}


def symbolic():
    import sympy as sp
    eta, slope = sp.symbols("eta slope", positive=True)
    return {"gamma": eta * slope}
