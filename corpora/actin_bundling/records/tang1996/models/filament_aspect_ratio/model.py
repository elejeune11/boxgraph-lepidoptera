def run(inputs, settings):
    return {"aspect_ratio": inputs["L"] / inputs["d_fil"]}


def symbolic():
    import sympy as sp

    L, d_fil = sp.symbols("L d_fil", positive=True)
    return {"aspect_ratio": L / d_fil}
