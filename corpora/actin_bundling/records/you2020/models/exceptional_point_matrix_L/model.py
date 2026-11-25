import math


def run(inputs, settings):
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    dc = math.sqrt(k * k + aB * aB)
    d = float(inputs["delta_over_dc"]) * dc
    kd_over_aB = (k * k - d * d) / aB
    # det[D_s] from Eq. (S24): (-kd/aB)*(-aB) - (-(k-d))*(-(k+d)) = kd - kd = 0
    det_Ds = (-kd_over_aB) * (-aB) - (-(k - d)) * (-(k + d))
    # Delta_1 = (aB + kd/aB)^2, so sqrt(Delta_1) = |aB + kd/aB|
    s_d1 = abs(aB + kd_over_aB)
    lam1_p = 0.5 * (-aB - kd_over_aB + s_d1)
    lam1_m = 0.5 * (-aB - kd_over_aB - s_d1)
    u1_p_1 = 0.5 * (aB - kd_over_aB + s_d1)
    u1_m_1 = 0.5 * (aB - kd_over_aB - s_d1)
    return {
        "det_Ds": det_Ds,
        "lambda_1_plus": lam1_p,
        "lambda_1_minus": lam1_m,
        "u1_plus_1": u1_p_1,
        "u1_minus_1": u1_m_1,
        "u1_plus_3": -u1_p_1,
        "u1_minus_3": -u1_m_1,
    }


def symbolic():
    import sympy as sp

    alpha_B, kappa, delta_over_dc = sp.symbols("alpha_B kappa delta_over_dc", positive=True)
    dc2 = kappa**2 + alpha_B**2
    d2 = delta_over_dc**2 * dc2
    kd_over_aB = (kappa**2 - d2) / alpha_B
    s_d1 = sp.Abs(alpha_B + kd_over_aB)
    u1_p_1 = sp.Rational(1, 2) * (alpha_B - kd_over_aB + s_d1)
    u1_m_1 = sp.Rational(1, 2) * (alpha_B - kd_over_aB - s_d1)
    return {
        "det_Ds": sp.Integer(0),
        "lambda_1_plus": sp.Rational(1, 2) * (-alpha_B - kd_over_aB + s_d1),
        "lambda_1_minus": sp.Rational(1, 2) * (-alpha_B - kd_over_aB - s_d1),
        "u1_plus_1": u1_p_1,
        "u1_minus_1": u1_m_1,
        "u1_plus_3": -u1_p_1,
        "u1_minus_3": -u1_m_1,
    }
