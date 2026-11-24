def run(inputs, settings):
    N_m = 1.0 / inputs["r_CP"]
    L_f = inputs["a_m"] * N_m
    return {"N_m": N_m, "L_f": L_f}


def symbolic():
    import sympy as sp

    a_m, r_CP = sp.symbols("a_m r_CP", positive=True)
    return {"N_m": 1 / r_CP, "L_f": a_m / r_CP}
