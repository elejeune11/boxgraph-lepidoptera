def run(inputs, settings):
    return {"tau_inv": 1.0 / (inputs["tau"] * 60.0)}


def symbolic():
    import sympy as sp
    tau = sp.symbols("tau", positive=True)
    return {"tau_inv": 1 / (60 * tau)}
