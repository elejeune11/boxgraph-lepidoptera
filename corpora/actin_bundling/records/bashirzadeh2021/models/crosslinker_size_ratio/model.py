def run(inputs, settings):
    return {"size_ratio": inputs["l_alpha"] / inputs["l_fascin"]}


def symbolic():
    import sympy as sp

    l_alpha, l_fascin = sp.symbols("l_alpha l_fascin", positive=True)
    return {"size_ratio": l_alpha / l_fascin}
