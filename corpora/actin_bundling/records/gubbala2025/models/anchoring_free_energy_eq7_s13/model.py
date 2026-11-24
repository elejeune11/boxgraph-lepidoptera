def run(inputs, settings):
    chi = float(inputs["chi"])
    S = float(inputs["S"])
    Q_xx = float(inputs["Q_xx"])
    Q_xy = float(inputs["Q_xy"])
    gx = float(inputs["dphi_x"])
    gy = float(inputs["dphi_y"])
    # Tensor term (Eq. 7): (chi / (2*S)) * Q_ij * d_i phi * d_j phi, with Q_yy = -Q_xx
    f_anch = (chi / (2.0 * S)) * (Q_xx * (gx * gx - gy * gy) + 2.0 * Q_xy * gx * gy)
    # Isotropic residual from Eq. (S13): (chi / 4) * |grad phi|^2
    f_iso = 0.25 * chi * (gx * gx + gy * gy)
    return {
        "f_anch": f_anch,
        "f_anch_iso_residual": f_iso,
        "f_anch_total": f_anch + f_iso,
    }


def symbolic():
    import sympy as sp

    chi, S, Q_xx, Q_xy, dphi_x, dphi_y = sp.symbols("chi S Q_xx Q_xy dphi_x dphi_y", real=True)
    f_anch = chi / (2 * S) * (Q_xx * (dphi_x**2 - dphi_y**2) + 2 * Q_xy * dphi_x * dphi_y)
    f_iso = chi / 4 * (dphi_x**2 + dphi_y**2)
    return {
        "f_anch": f_anch,
        "f_anch_iso_residual": f_iso,
        "f_anch_total": sp.simplify(f_anch + f_iso),
    }
