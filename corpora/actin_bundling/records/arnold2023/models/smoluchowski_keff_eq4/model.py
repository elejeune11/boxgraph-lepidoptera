import math

K_B = 1.380649e-23  # J / K


def run(inputs, settings):
    T = float(inputs["T"])
    eta = float(inputs["eta"])
    a_um = float(inputs["a"])
    L_um = float(inputs["L"])
    gamma_dot = float(inputs["gamma_dot"])
    C = float(inputs["C"])
    n = float(inputs["n"])

    a_m = a_um * 1.0e-6
    D_c_m2s = K_B * T / (16.0 * eta * a_m)
    Pe_c = gamma_dot * (a_m ** 2) / D_c_m2s

    k_eff_m2s = (math.pi * K_B * T / (8.0 * eta * a_m * math.log(L_um / (2.0 * a_um)))) * (
        1.0 + C * math.sqrt(Pe_c)
    )
    k_eff_um2s = k_eff_m2s * 1.0e12
    J = k_eff_um2s * (n ** 2)

    # Eq. (22) / Eq. (32) normalized pair distribution at geometric midpoint r = sqrt(2*a*L)
    r_mid = math.sqrt(2.0 * a_um * L_um)
    g_norm_mid = math.log(r_mid / (2.0 * a_um)) / math.log(L_um / (2.0 * a_um))
    return {
        "k_eff": float(k_eff_um2s),
        "J": float(J),
        "Pe_c": float(Pe_c),
        "g_norm_mid": float(g_norm_mid),
    }


def symbolic():
    import sympy as sp
    T, eta, a, L, gamma_dot, C, n = sp.symbols("T eta a L gamma_dot C n", positive=True)
    D_c_um2s = K_B * T * 1e18 / (16 * eta * a)
    Pe_c = gamma_dot * a**2 / D_c_um2s
    k_eff = (sp.pi * K_B * T * 1e18 / (8 * eta * a * sp.log(L / (2 * a)))) * (1 + C * sp.sqrt(Pe_c))
    J = k_eff * n**2
    g_norm_mid = sp.Rational(1, 2)
    return {"k_eff": k_eff, "J": J, "Pe_c": Pe_c, "g_norm_mid": g_norm_mid}
