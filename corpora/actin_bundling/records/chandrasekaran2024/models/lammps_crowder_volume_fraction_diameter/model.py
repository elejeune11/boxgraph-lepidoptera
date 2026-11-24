import math


def run(inputs, settings):
    v_sphere = (math.pi / 6.0) * (inputs["d_crowder"] ** 3)
    v_box = inputs["L_box"] ** 3
    return {"phi": 100.0 * inputs["N_crowder"] * v_sphere / v_box}


def symbolic():
    import sympy as sp

    N_crowder, d_crowder, L_box = sp.symbols("N_crowder d_crowder L_box", positive=True)
    return {"phi": 100 * N_crowder * (sp.pi / 6) * d_crowder**3 / L_box**3}
