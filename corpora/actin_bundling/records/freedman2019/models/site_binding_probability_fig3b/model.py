import math


def run(inputs, settings):
    n_g = inputs["n_g"]
    l_b = inputs["l_b"]
    L_p_nm = inputs["L_p"] * 1000.0
    dl = inputs["delta_l_xl"]
    l_g_nm = n_g * l_b
    l_g_um = l_g_nm * 1e-3
    if l_g_nm <= dl:
        return {"l_g": l_g_um, "U_kBT": float("inf"), "P_diff": 0.0}
    theta = math.asin(dl / l_g_nm)
    u_kbt = (L_p_nm / l_g_nm) * (theta * theta)
    if u_kbt > 700.0:
        p_diff = 0.0
    else:
        e = math.exp(-u_kbt)
        p_diff = e / (1.0 + e)
    return {"l_g": l_g_um, "U_kBT": u_kbt, "P_diff": p_diff}


def symbolic():
    import sympy as sp

    n_g, l_b, L_p, delta_l_xl = sp.symbols("n_g l_b L_p delta_l_xl", positive=True)
    l_g_nm = n_g * l_b
    u_kbt = (1000 * L_p / l_g_nm) * sp.asin(delta_l_xl / l_g_nm) ** 2
    e = sp.exp(-u_kbt)
    return {
        "l_g": l_g_nm * sp.Rational(1, 1000),
        "U_kBT": u_kbt,
        "P_diff": e / (1 + e),
    }
