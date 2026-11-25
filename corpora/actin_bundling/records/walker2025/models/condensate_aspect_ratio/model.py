def run(inputs, settings):
    return {"AR": float(inputs["L_long"]) / float(inputs["L_short"])}


def symbolic():
    import sympy as sp

    L_long, L_short = sp.symbols("L_long L_short", positive=True)
    return {"AR": L_long / L_short}
