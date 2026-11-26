def run(inputs, settings):
    return {"c_lys18_total": inputs["c_18"] + inputs["s"] * inputs["c_actin"]}


def symbolic():
    import sympy as sp

    c_18, s, c_actin = sp.symbols("c_18 s c_actin", positive=True)
    return {"c_lys18_total": c_18 + s * c_actin}
