def run(inputs, settings):
    return {"RAM_crit": inputs["c_star"] / inputs["c_global"]}


def symbolic():
    import sympy as sp

    c_star, c_global = sp.symbols("c_star c_global", positive=True)
    return {"RAM_crit": c_star / c_global}
