import math

NA_UM3_PER_UMOL = 602.214076


def run(inputs, settings):
    r_um = float(inputs["r_bind"]) * 1e-3
    L_um = float(inputs["L_seg"])
    V_bind = math.pi * r_um * r_um * L_um
    k_bind_2 = float(inputs["k_bind"]) * NA_UM3_PER_UMOL * V_bind
    return {"k_bind_2": k_bind_2, "V_bind": V_bind}


def symbolic():
    import sympy as sp

    k_bind, r_bind, L_seg = sp.symbols("k_bind r_bind L_seg", positive=True)
    r_um = r_bind * sp.Rational(1, 1000)
    V_bind = sp.pi * r_um**2 * L_seg
    k_bind_2 = k_bind * sp.Float(NA_UM3_PER_UMOL) * V_bind
    return {"k_bind_2": k_bind_2, "V_bind": V_bind}
