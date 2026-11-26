def run(inputs, settings):
    theta = inputs["c_ABP"] / (inputs["K_d"] + inputs["c_ABP"])
    delta = inputs["delta_max"] / theta
    return {"theta": theta, "delta": delta}


def symbolic():
    import sympy as sp

    c_ABP, K_d, delta_max = sp.symbols("c_ABP K_d delta_max", positive=True)
    theta = c_ABP / (K_d + c_ABP)
    return {"theta": theta, "delta": delta_max / theta}
