def run(inputs, settings):
    return {"omega": (inputs["aspect_ratio"] ** 2) / 4.0}


def symbolic():
    import sympy as sp

    aspect_ratio = sp.Symbol("aspect_ratio", positive=True)
    return {"omega": aspect_ratio**2 / 4}
