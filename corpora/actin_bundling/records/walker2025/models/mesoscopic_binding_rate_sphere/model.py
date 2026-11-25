import math

NA_UM3_PER_UMOL = 602.214076


def run(inputs, settings):
    r_um = float(inputs["r_bind"]) * 1e-3
    V_bind = (4.0 / 3.0) * math.pi * (r_um ** 3)
    k_bind = float(inputs["k_bind_2"]) / (NA_UM3_PER_UMOL * V_bind)
    return {"k_bind": k_bind, "V_bind": V_bind}


def symbolic():
    import sympy as sp

    k_bind_2, r_bind = sp.symbols("k_bind_2 r_bind", positive=True)
    r_um = r_bind * sp.Rational(1, 1000)
    V_bind = sp.Rational(4, 3) * sp.pi * r_um**3
    k_bind = k_bind_2 / (sp.Float(NA_UM3_PER_UMOL) * V_bind)
    return {"k_bind": k_bind, "V_bind": V_bind}
