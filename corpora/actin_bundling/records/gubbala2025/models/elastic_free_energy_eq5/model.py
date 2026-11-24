def run(inputs, settings):
    E1 = float(inputs["E1"])
    E3 = float(inputs["E3"])
    Q_xx = float(inputs["Q_xx"])
    Q_xy = float(inputs["Q_xy"])
    qx_x = float(inputs["dQxx_x"])
    qx_y = float(inputs["dQxx_y"])
    qxy_x = float(inputs["dQxy_x"])
    qxy_y = float(inputs["dQxy_y"])
    # For 2D symmetric traceless Q: (d_i Q_kl)(d_j Q_kl) = 2*(d_i Q_xx d_j Q_xx + d_i Q_xy d_j Q_xy)
    g_xx = 2.0 * (qx_x * qx_x + qxy_x * qxy_x)
    g_yy = 2.0 * (qx_y * qx_y + qxy_y * qxy_y)
    g_xy = 2.0 * (qx_x * qx_y + qxy_x * qxy_y)
    f_E1 = 0.5 * E1 * (g_xx + g_yy)
    f_E3 = 0.5 * E3 * (Q_xx * (g_xx - g_yy) + 2.0 * Q_xy * g_xy)
    return {"f_elastic": f_E1 + f_E3, "f_E1": f_E1, "f_E3": f_E3}


def symbolic():
    import sympy as sp

    E1, E3, Q_xx, Q_xy, dQxx_x, dQxx_y, dQxy_x, dQxy_y = sp.symbols(
        "E1 E3 Q_xx Q_xy dQxx_x dQxx_y dQxy_x dQxy_y", real=True
    )
    g_xx = 2 * (dQxx_x**2 + dQxy_x**2)
    g_yy = 2 * (dQxx_y**2 + dQxy_y**2)
    g_xy = 2 * (dQxx_x * dQxx_y + dQxy_x * dQxy_y)
    f_E1 = sp.Rational(1, 2) * E1 * (g_xx + g_yy)
    f_E3 = sp.Rational(1, 2) * E3 * (Q_xx * (g_xx - g_yy) + 2 * Q_xy * g_xy)
    return {"f_elastic": sp.expand(f_E1 + f_E3), "f_E1": sp.expand(f_E1), "f_E3": sp.expand(f_E3)}
