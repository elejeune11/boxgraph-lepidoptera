def run(inputs, settings):
    return {"P_disk": 100.0 - inputs["P_ring"] - inputs["P_shell"]}


def symbolic():
    import sympy as sp

    P_ring, P_shell = sp.symbols("P_ring P_shell", positive=True)
    return {"P_disk": 100 - P_ring - P_shell}
