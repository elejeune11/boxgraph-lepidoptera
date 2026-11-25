import math

NA_UM3_PER_UMOL = 602.214076  # N_A [1/mol] * 1e-15 [L/um^3] * 1e-6 [mol/umol]


def run(inputs, settings):
    r_um = float(inputs["r_bind"]) * 1e-3
    L_um = float(inputs["L_seg"])
    V_bind = math.pi * r_um * r_um * L_um  # um^3
    k_bind = float(inputs["k_bind_2"]) / (NA_UM3_PER_UMOL * V_bind)
    return {"k_bind": k_bind, "V_bind": V_bind}


def symbolic():
    import sympy as sp

    k_bind_2, r_bind, L_seg = sp.symbols("k_bind_2 r_bind L_seg", positive=True)
    r_um = r_bind * sp.Rational(1, 1000)
    V_bind = sp.pi * r_um**2 * L_seg
    k_bind = k_bind_2 / (sp.Float(NA_UM3_PER_UMOL) * V_bind)
    return {"k_bind": k_bind, "V_bind": V_bind}
