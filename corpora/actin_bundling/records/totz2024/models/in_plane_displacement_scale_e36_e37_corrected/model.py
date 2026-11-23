def run(inputs, settings):
    L = inputs["L"]
    h = inputs["h"]
    return {"u_scale": (h ** 2) / (L / 2.0)}


def symbolic():
    import sympy as sp
    L, h = sp.symbols("L h", positive=True)
    return {"u_scale": h**2 / (L / 2)}
