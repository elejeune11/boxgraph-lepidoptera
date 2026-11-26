def run(inputs, settings):
    I_t = float(inputs["I_t"])
    I_inf = float(inputs["I_inf"])
    return {"c_0": 1.0 - I_t / I_inf}


def symbolic():
    import sympy as sp

    I_t, I_inf = sp.symbols("I_t I_inf", positive=True)
    return {"c_0": 1 - I_t / I_inf}
