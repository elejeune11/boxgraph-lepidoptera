def run(inputs, settings):
    return {"ratio": inputs["t0_fb5"] / inputs["tau_c_fb2"]}


def symbolic():
    import sympy as sp

    t0_fb5, tau_c_fb2 = sp.symbols("t0_fb5 tau_c_fb2", positive=True)
    return {"ratio": t0_fb5 / tau_c_fb2}
