def run(inputs, settings):
    k = inputs["k"]          # pN/um
    Gamma = inputs["Gamma"]  # pN*s/um
    d = inputs["d"]          # um
    tau_s = Gamma / k        # s
    v_um_s = (k / Gamma) * d
    return {
        "tau": tau_s / 60.0,        # min
        "v": v_um_s * 1000.0,       # nm/s
    }


def symbolic():
    import sympy as sp
    k, Gamma, d = sp.symbols("k Gamma d", positive=True)
    return {
        "tau": Gamma / (60 * k),
        "v": 1000 * k * d / Gamma,
    }
