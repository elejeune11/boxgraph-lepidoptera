def run(inputs, settings):
    return {"R_mob": inputs["t_half_ml"] / inputs["t_half_slb"]}


def symbolic():
    import sympy as sp

    t_half_ml, t_half_slb = sp.symbols("t_half_ml t_half_slb", positive=True)
    return {"R_mob": t_half_ml / t_half_slb}
