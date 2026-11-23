def run(inputs, settings):
    return {"R_star": 1.0 / inputs["r_sat"]}


def symbolic():
    import sympy as sp

    r_sat = sp.Symbol("r_sat", positive=True)
    return {"R_star": 1 / r_sat}
