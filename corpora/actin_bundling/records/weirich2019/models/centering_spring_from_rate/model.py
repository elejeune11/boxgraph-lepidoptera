def run(inputs, settings):
    return {"k": inputs["Gamma"] * inputs["tau_inv"]}


def symbolic():
    import sympy as sp
    Gamma, tau_inv = sp.symbols("Gamma tau_inv", positive=True)
    return {"k": Gamma * tau_inv}
