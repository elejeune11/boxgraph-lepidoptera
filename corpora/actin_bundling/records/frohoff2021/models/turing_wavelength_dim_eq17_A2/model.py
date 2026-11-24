import math


def run(inputs, settings):
    L0 = inputs["L0"]
    f2_pp = inputs["f2_pp"]
    Q = inputs["Q"]
    kappa = inputs["kappa"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    Delta = alpha * alpha - rho * rho
    if Delta <= 0.0 or kappa <= 0.0 or L0 <= 0.0 or abs(Q * kappa - 1.0) < 1e-15:
        return {
            "qc_dim": None,
            "inv_qc_dim": None,
            "lambda_c_dim": None,
            "lambda_c_dimless": None,
        }

    sign = 1.0 if kappa > 1.0 / Q else -1.0
    qc_sq = (sign * math.sqrt(kappa * Delta) - f2_pp) / kappa
    if qc_sq <= 0.0:
        return {
            "qc_dim": None,
            "inv_qc_dim": None,
            "lambda_c_dim": None,
            "lambda_c_dimless": None,
        }

    qc = math.sqrt(qc_sq)
    return {
        "qc_dim": qc / L0,
        "inv_qc_dim": L0 / qc,
        "lambda_c_dim": 2.0 * math.pi * L0 / qc,
        "lambda_c_dimless": 2.0 * math.pi / qc,
    }


def symbolic():
    import sympy as sp

    L0, f2_pp, Q, kappa, rho, alpha = sp.symbols("L0 f2_pp Q kappa rho alpha", real=True)
    Delta = alpha**2 - rho**2
    sign = sp.Piecewise((1, kappa * Q > 1), (-1, True))
    qc = sp.sqrt((sign * sp.sqrt(kappa * Delta) - f2_pp) / kappa)
    return {
        "qc_dim": qc / L0,
        "inv_qc_dim": L0 / qc,
        "lambda_c_dim": 2 * sp.pi * L0 / qc,
        "lambda_c_dimless": 2 * sp.pi / qc,
    }
