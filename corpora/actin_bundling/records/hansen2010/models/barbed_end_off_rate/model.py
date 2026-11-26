def run(inputs, settings):
    return {"K_off": 1.0 / inputs["tau_1"]}


def symbolic():
    import sympy as sp

    tau_1 = sp.symbols("tau_1", positive=True)
    return {"K_off": 1 / tau_1}
