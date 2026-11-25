import cmath
import math


def run(inputs, settings):
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    chi = float(inputs["chi"])
    alpha = float(inputs["alpha"])
    q = float(inputs["q"])
    q2 = q * q

    D11 = 4.0 * q2 * c1 * c1
    D12 = -q2 * (chi + alpha)
    D21 = -q2 * (chi - alpha)
    D22 = 4.0 * q2 * c2 * c2

    alpha_star = math.sqrt(chi * chi + 4.0 * (c1 * c1 - c2 * c2) ** 2)
    rad = (alpha_star + alpha) * (alpha_star - alpha)
    if rad >= 0.0:
        sr = math.sqrt(rad)
        lam1 = 2.0 * q2 * (c1 * c1 + c2 * c2) + q2 * sr
        lam2 = 2.0 * q2 * (c1 * c1 + c2 * c2) - q2 * sr
        re_max = lam1
        im_abs = 0.0
    else:
        si = math.sqrt(-rad)
        lam1 = complex(2.0 * q2 * (c1 * c1 + c2 * c2), q2 * si)
        lam2 = complex(2.0 * q2 * (c1 * c1 + c2 * c2), -q2 * si)
        re_max = lam1.real
        im_abs = q2 * si

    # Eigenvectors from Eq. (17): eta_k = [lam_k - D22, D21]^T
    v1_0, v1_1 = complex(lam1 - D22), complex(D21)
    v2_0, v2_1 = complex(lam2 - D22), complex(D21)
    n1 = math.sqrt(abs(v1_0) ** 2 + abs(v1_1) ** 2)
    n2 = math.sqrt(abs(v2_0) ** 2 + abs(v2_1) ** 2)
    if n1 > 0.0 and n2 > 0.0:
        cross = abs(v1_0 * v2_1 - v1_1 * v2_0)
        sin_ang = min(1.0, cross / (n1 * n2))
    else:
        sin_ang = 0.0

    return {
        "D11": D11,
        "D12": D12,
        "D21": D21,
        "D22": D22,
        "alpha_star": alpha_star,
        "re_lambda_max": re_max,
        "im_lambda_abs": im_abs,
        "eigvec_sin_angle": sin_ang,
    }


def symbolic():
    import sympy as sp

    c1, c2, chi, alpha, q = sp.symbols("c1 c2 chi alpha q", real=True)
    return {
        "D11": 4 * q**2 * c1**2,
        "D12": -q**2 * (chi + alpha),
        "D21": -q**2 * (chi - alpha),
        "D22": 4 * q**2 * c2**2,
        "alpha_star": sp.sqrt(chi**2 + 4 * (c1**2 - c2**2) ** 2),
    }
