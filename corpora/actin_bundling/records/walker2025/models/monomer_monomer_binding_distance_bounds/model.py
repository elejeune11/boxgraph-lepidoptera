def run(inputs, settings):
    R = float(inputs["R_solid"])
    return {"d_min": R, "d_max": 3.0 * R}


def symbolic():
    import sympy as sp

    R_solid = sp.Symbol("R_solid", positive=True)
    return {"d_min": R_solid, "d_max": 3 * R_solid}
