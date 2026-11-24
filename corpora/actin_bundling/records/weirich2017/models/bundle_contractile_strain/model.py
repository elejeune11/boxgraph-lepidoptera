def run(inputs, settings):
    return {"strain": (inputs["L_i"] - inputs["L_f"]) / inputs["L_i"]}


def symbolic():
    import sympy as sp

    L_i, L_f = sp.symbols("L_i L_f", positive=True)
    return {"strain": (L_i - L_f) / L_i}
