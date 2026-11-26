import math


def run(inputs, settings):
    k1 = float(inputs["k1"])
    k2 = float(inputs["k2"])
    t = float(inputs["t"])
    p_atp = math.exp(-k1 * t)
    if abs(k1 - k2) < 1e-12:
        p_adp_pi = k1 * t * math.exp(-k1 * t)
    else:
        p_adp_pi = (k1 / (k1 - k2)) * (math.exp(-k2 * t) - math.exp(-k1 * t))
    p_adp = 1.0 - p_atp - p_adp_pi
    return {"p_atp": p_atp, "p_adp_pi": p_adp_pi, "p_adp": p_adp}


def symbolic():
    import sympy as sp

    k1, k2, t = sp.symbols("k1 k2 t", positive=True)
    p_atp = sp.exp(-k1 * t)
    p_adp_pi = (k1 / (k1 - k2)) * (sp.exp(-k2 * t) - sp.exp(-k1 * t))
    p_adp = 1 - p_atp - p_adp_pi
    return {"p_atp": p_atp, "p_adp_pi": p_adp_pi, "p_adp": p_adp}
