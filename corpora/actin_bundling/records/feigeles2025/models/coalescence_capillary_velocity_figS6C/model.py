def run(inputs, settings):
    # (1 uN/m) / (1 Pa*s) = (1e-6 N/m) / (N*s/m^2) = 1e-6 m/s = 1 um/s
    # (1 Pa*s) / (1 uN/m) = 1 s/um
    eta = inputs["eta"]
    gamma = inputs["gamma"]
    return {
        "v_cap": gamma / eta,
        "inv_v_cap": eta / gamma,
    }


def symbolic():
    import sympy as sp
    eta, gamma = sp.symbols("eta gamma", positive=True)
    return {
        "v_cap": gamma / eta,
        "inv_v_cap": eta / gamma,
    }
