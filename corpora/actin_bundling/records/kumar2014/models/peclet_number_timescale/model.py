def run(inputs, settings):
    return {"P": inputs["tau"] / inputs["tau_a"]}


def symbolic():
    import sympy as sp

    tau, tau_a = sp.symbols("tau tau_a", positive=True)
    return {"P": tau / tau_a}
