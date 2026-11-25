def run(inputs, settings):
    return {"L": (float(inputs["N_x"]) - 1.0) * float(inputs["h"])}


def symbolic():
    import sympy as sp

    N_x, h = sp.symbols("N_x h", positive=True)
    return {"L": (N_x - 1) * h}
