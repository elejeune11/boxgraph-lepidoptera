def run(inputs, settings):
    return {"r_mon": inputs["v"] / inputs["delta_l"]}


def symbolic():
    import sympy as sp

    v, delta_l = sp.symbols("v delta_l", positive=True)
    return {"r_mon": v / delta_l}
