def run(inputs, settings):
    # (1 ms / 1 um) * (1 mN / 1 m) = (1e-3 s / 1e-6 m) * (1e-3 N / m) = 1 Pa*s
    return {"eta": inputs["v_inv"] * inputs["gamma"]}


def symbolic():
    import sympy as sp

    v_inv, gamma = sp.symbols("v_inv gamma", positive=True)
    return {"eta": v_inv * gamma}
