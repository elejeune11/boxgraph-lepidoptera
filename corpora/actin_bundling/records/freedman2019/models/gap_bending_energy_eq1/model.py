import math

# Boltzmann constant in pN*um/K: 1.380649e-23 J/K * 1e18 (pN*um)/J
KB_PN_UM_PER_K = 1.380649e-5


def run(inputs, settings):
    L_p_nm = inputs["L_p"] * 1000.0
    dl_nm = inputs["delta_l_xl"]
    l_g_nm = inputs["l_g"]
    T = inputs["T"]
    theta = math.asin(dl_nm / l_g_nm)
    u_kbt = (L_p_nm / l_g_nm) * (theta * theta)
    return {
        "U_kBT": u_kbt,
        "U": u_kbt * KB_PN_UM_PER_K * T,
    }


def symbolic():
    import sympy as sp

    L_p, delta_l_xl, l_g, T = sp.symbols("L_p delta_l_xl l_g T", positive=True)
    u_kbt = (1000 * L_p / l_g) * sp.asin(delta_l_xl / l_g) ** 2
    return {
        "U_kBT": u_kbt,
        "U": u_kbt * KB_PN_UM_PER_K * T,
    }
