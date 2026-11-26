import math


def run(inputs, settings):
    l_g_nm = inputs["n_g"] * inputs["l_b"]
    dl = inputs["delta_l_xl"]
    L_p_nm = inputs["L_p"] * 1000.0
    u_gap = (L_p_nm / l_g_nm) * (math.asin(dl / l_g_nm) ** 2)
    u_lat = -(inputs["N_S"] * inputs["mu_S_kBT"] + inputs["N_L"] * inputs["mu_L_kBT"]) + u_gap
    return {"U_lat_kBT": u_lat}


def symbolic():
    import sympy as sp

    N_S, N_L, n_g, l_b, L_p, delta_l_xl = sp.symbols(
        "N_S N_L n_g l_b L_p delta_l_xl", positive=True
    )
    mu_S_kBT, mu_L_kBT = sp.symbols("mu_S_kBT mu_L_kBT", real=True)
    l_g_nm = n_g * l_b
    u_gap = (1000 * L_p / l_g_nm) * sp.asin(delta_l_xl / l_g_nm) ** 2
    return {"U_lat_kBT": -(N_S * mu_S_kBT + N_L * mu_L_kBT) + u_gap}
