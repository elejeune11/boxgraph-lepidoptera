def run(inputs, settings):
    return {"alpha_star": abs(float(inputs["chi"]))}


def symbolic():
    import sympy as sp

    chi = sp.Symbol("chi", real=True)
    return {"alpha_star": sp.Abs(chi)}
