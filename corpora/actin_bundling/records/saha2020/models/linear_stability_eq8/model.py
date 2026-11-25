import math


def run(inputs, settings):
    p1 = float(inputs["phi1_bar"])
    p2 = float(inputs["phi2_bar"])
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    alpha = float(inputs["alpha"])
    q = float(inputs["q"])
    q2 = q * q

    d11 = -(2.0 * (p1 - c11) ** 2 + 8.0 * (p1 - c11) * (p1 - c12) + 2.0 * (p1 - c12) ** 2 + 2.0 * p2 * p2 * chip)
    d22 = -(2.0 * (p2 - c21) ** 2 + 8.0 * (p2 - c21) * (p2 - c22) + 2.0 * (p2 - c22) ** 2 + 2.0 * p1 * p1 * chip)
    d12 = -((chi + alpha) + 4.0 * p1 * p2 * chip)
    d21 = -((chi - alpha) + 4.0 * p1 * p2 * chip)

    tr_over_q2 = d11 + d22
    disc_over_q4 = (d11 - d22) ** 2 + 4.0 * d12 * d21
    if disc_over_q4 >= 0.0:
        re_max = 0.5 * q2 * (tr_over_q2 + math.sqrt(disc_over_q4))
        im_abs = 0.0
    else:
        re_max = 0.5 * q2 * tr_over_q2
        im_abs = 0.5 * q2 * math.sqrt(-disc_over_q4)

    return {
        "D11": q2 * d11,
        "D12": q2 * d12,
        "D21": q2 * d21,
        "D22": q2 * d22,
        "trace_D": q2 * tr_over_q2,
        "disc_over_q4": disc_over_q4,
        "re_lambda_max": re_max,
        "im_lambda_abs": im_abs,
    }


def symbolic():
    import sympy as sp

    phi1_bar, phi2_bar, c11, c12, c21, c22, chi, chip, alpha, q = sp.symbols(
        "phi1_bar phi2_bar c11 c12 c21 c22 chi chip alpha q", real=True
    )
    d11 = -(2 * (phi1_bar - c11) ** 2 + 8 * (phi1_bar - c11) * (phi1_bar - c12) + 2 * (phi1_bar - c12) ** 2 + 2 * phi2_bar**2 * chip)
    d22 = -(2 * (phi2_bar - c21) ** 2 + 8 * (phi2_bar - c21) * (phi2_bar - c22) + 2 * (phi2_bar - c22) ** 2 + 2 * phi1_bar**2 * chip)
    d12 = -((chi + alpha) + 4 * phi1_bar * phi2_bar * chip)
    d21 = -((chi - alpha) + 4 * phi1_bar * phi2_bar * chip)
    disc_over_q4 = (d11 - d22) ** 2 + 4 * d12 * d21
    return {
        "D11": q**2 * d11,
        "D12": q**2 * d12,
        "D21": q**2 * d21,
        "D22": q**2 * d22,
        "trace_D": q**2 * (d11 + d22),
        "disc_over_q4": disc_over_q4,
    }
