def run(inputs, settings):
    chi = float(inputs["chi"])
    return {"lam": 1.0 + 0.5 * abs(chi)}


def symbolic():
    import sympy as sp

    chi = sp.symbols("chi", real=True)
    return {"lam": 1 + sp.Abs(chi) / 2}
