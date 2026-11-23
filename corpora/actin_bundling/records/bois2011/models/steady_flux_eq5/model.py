def run(inputs, settings):
    # From Eq. (5), v = D * d_x(ln c) + j0 / c = (D / c) * d_x c + j0 / c,
    # so j0 = c * v - D * d_x c.
    return {"j0": inputs["c"] * inputs["v"] - inputs["D"] * inputs["dc_dx"]}


def symbolic():
    import sympy as sp

    v, dc_dx = sp.symbols("v dc_dx", real=True)
    c, D = sp.symbols("c D", positive=True)
    return {"j0": c * v - D * dc_dx}
