import math


def run(inputs, settings):
    v_sphere = (4.0 / 3.0) * math.pi * (inputs["r_crowder"] ** 3)
    v_box = inputs["L_box"] ** 3
    return {"phi": 100.0 * inputs["N_crowder"] * v_sphere / v_box}


def symbolic():
    import sympy as sp

    N_crowder, r_crowder, L_box = sp.symbols("N_crowder r_crowder L_box", positive=True)
    return {"phi": 100 * N_crowder * sp.Rational(4, 3) * sp.pi * r_crowder**3 / L_box**3}
