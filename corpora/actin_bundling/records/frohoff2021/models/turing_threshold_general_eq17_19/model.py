import math


def run(inputs, settings):
    f2_pp = inputs["f2_pp"]
    Q = inputs["Q"]
    kappa = inputs["kappa"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]

    Delta = alpha * alpha - rho * rho
    if Delta <= 0.0 or kappa <= 0.0 or abs(Q * kappa - 1.0) < 1e-15:
        return {"f1_pp_T": None, "qc_sq": None, "qc": None, "f2_pp_T_end": None}

    sign = 1.0 if kappa > 1.0 / Q else -1.0
    skd = math.sqrt(kappa * Delta)
    f2_end = sign * skd
    qc_sq = (f2_end - f2_pp) / kappa
    if qc_sq <= 0.0:
        return {
            "f1_pp_T": (f2_pp - 2.0 * f2_end) / kappa,
            "qc_sq": qc_sq,
            "qc": None,
            "f2_pp_T_end": f2_end,
        }
    return {
        "f1_pp_T": (f2_pp - 2.0 * f2_end) / kappa,
        "qc_sq": qc_sq,
        "qc": math.sqrt(qc_sq),
        "f2_pp_T_end": f2_end,
    }


def symbolic():
    import sympy as sp

    f2_pp, Q, kappa, rho, alpha = sp.symbols("f2_pp Q kappa rho alpha", real=True)
    Delta = alpha**2 - rho**2
    skd = sp.sqrt(kappa * Delta)
    sign = sp.Piecewise((1, kappa * Q > 1), (-1, True))
    f2_end = sign * skd
    qc_sq = (f2_end - f2_pp) / kappa
    return {
        "f1_pp_T": (f2_pp - 2 * f2_end) / kappa,
        "qc_sq": qc_sq,
        "qc": sp.sqrt(qc_sq),
        "f2_pp_T_end": f2_end,
    }
