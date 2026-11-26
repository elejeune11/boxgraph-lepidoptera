def run(inputs, settings):
    return {"N": inputs["v"] * inputs["tau_1"]}


def symbolic():
    import sympy as sp

    v, tau_1 = sp.symbols("v tau_1", positive=True)
    return {"N": v * tau_1}
