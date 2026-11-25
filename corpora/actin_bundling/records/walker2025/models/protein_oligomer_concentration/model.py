def run(inputs, settings):
    return {"C_olig": float(inputs["C_mon"]) / float(inputs["n_subunits"])}


def symbolic():
    import sympy as sp

    C_mon, n_subunits = sp.symbols("C_mon n_subunits", positive=True)
    return {"C_olig": C_mon / n_subunits}
