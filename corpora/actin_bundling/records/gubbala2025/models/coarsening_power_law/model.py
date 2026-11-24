def run(inputs, settings):
    t = float(inputs["t"])
    t0 = float(inputs["t0"])
    alpha = float(inputs["alpha"])
    return {"R_over_R0": (t / t0) ** alpha}


def symbolic():
    import sympy as sp

    t, t0 = sp.symbols("t t0", positive=True)
    alpha = sp.symbols("alpha", real=True)
    return {"R_over_R0": (t / t0) ** alpha}
