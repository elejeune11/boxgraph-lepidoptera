def run(inputs, settings):
    return {"phi": (inputs["l_h"] / inputs["l_t"]) ** 2}


def symbolic():
    import sympy as sp

    l_h, l_t = sp.symbols("l_h l_t", positive=True)
    return {"phi": (l_h / l_t) ** 2}
