def run(inputs, settings):
    return {"T_star": inputs["kBT"] / inputs["eps_VV"]}


def symbolic():
    import sympy as sp

    kBT, eps_VV = sp.symbols("kBT eps_VV", positive=True)
    return {"T_star": kBT / eps_VV}
