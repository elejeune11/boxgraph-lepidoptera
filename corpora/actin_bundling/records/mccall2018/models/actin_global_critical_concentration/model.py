def run(inputs, settings):
    return {"c_global_crit": inputs["c_star"] / inputs["RAM"]}


def symbolic():
    import sympy as sp

    c_star, RAM = sp.symbols("c_star RAM", positive=True)
    return {"c_global_crit": c_star / RAM}
