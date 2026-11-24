def run(inputs, settings):
    return {"N_fil": inputs["T_ring"] / inputs["delta"]}


def symbolic():
    import sympy as sp

    T_ring, delta = sp.symbols("T_ring delta", positive=True)
    return {"N_fil": T_ring / delta}
