def run(inputs, settings):
    # k_grow in nm/s -> 1e-3 um/s
    return {"L_fil": inputs["L0"] + inputs["k_grow"] * 1e-3 * inputs["t"]}


def symbolic():
    import sympy as sp

    L0, k_grow, t = sp.symbols("L0 k_grow t", positive=True)
    return {"L_fil": L0 + sp.Rational(1, 1000) * k_grow * t}
