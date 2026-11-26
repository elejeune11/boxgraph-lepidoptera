def run(inputs, settings):
    return {"delta_l_xl": inputs["l_long"] - inputs["l_short"]}


def symbolic():
    import sympy as sp

    l_long, l_short = sp.symbols("l_long l_short", positive=True)
    return {"delta_l_xl": l_long - l_short}
