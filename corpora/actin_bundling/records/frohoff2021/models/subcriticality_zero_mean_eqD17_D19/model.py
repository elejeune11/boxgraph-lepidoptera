import math


def run(inputs, settings):
    q_n = inputs["q_n"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    eps = inputs["eps"]

    M_n = q_n * q_n * (1.0 - kappa) - a_Delta
    Delta = alpha * alpha - rho * rho
    Sigma = (rho + alpha) ** 2
    Delta_H_n = 0.25 * M_n * M_n
    alpha_H_n = math.sqrt(rho * rho + Delta_H_n)

    rad = Delta_H_n - Delta
    if rad < 0.0:
        return {
            "M_n": M_n,
            "xi_plus": None,
            "xi_minus": None,
            "subcrit_ratio_plus": None,
            "subcrit_ratio_minus": None,
            "Delta_H_n": Delta_H_n,
            "alpha_H_n": alpha_H_n,
            "A0_plus": None,
        }

    s = math.sqrt(rad)
    xi_plus = (0.5 * M_n + s) ** 2
    xi_minus = (0.5 * M_n - s) ** 2

    den_p = -Delta * Sigma + xi_plus * xi_plus
    den_m = -Delta * Sigma + xi_minus * xi_minus
    R_p = (Delta - xi_plus) / den_p if abs(den_p) > 1e-15 else None
    R_m = (Delta - xi_minus) / den_m if abs(den_m) > 1e-15 else None

    A0_p = None
    if R_p is not None:
        val = (eps / (3.0 * math.sqrt(2.0))) * R_p
        if val >= 0.0:
            A0_p = math.sqrt(val)

    return {
        "M_n": M_n,
        "xi_plus": xi_plus,
        "xi_minus": xi_minus,
        "subcrit_ratio_plus": R_p,
        "subcrit_ratio_minus": R_m,
        "Delta_H_n": Delta_H_n,
        "alpha_H_n": alpha_H_n,
        "A0_plus": A0_p,
    }


def symbolic():
    import sympy as sp

    q_n, a_Delta, kappa, rho, alpha, eps = sp.symbols(
        "q_n a_Delta kappa rho alpha eps", real=True
    )
    M_n = q_n**2 * (1 - kappa) - a_Delta
    Delta = alpha**2 - rho**2
    Sigma = (rho + alpha) ** 2
    Delta_H_n = M_n**2 / 4
    s = sp.sqrt(Delta_H_n - Delta)
    xi_plus = (M_n / 2 + s) ** 2
    xi_minus = (M_n / 2 - s) ** 2
    R_p = (Delta - xi_plus) / (-Delta * Sigma + xi_plus**2)
    R_m = (Delta - xi_minus) / (-Delta * Sigma + xi_minus**2)
    return {
        "M_n": M_n,
        "xi_plus": xi_plus,
        "xi_minus": xi_minus,
        "subcrit_ratio_plus": R_p,
        "subcrit_ratio_minus": R_m,
        "Delta_H_n": Delta_H_n,
        "alpha_H_n": sp.sqrt(rho**2 + Delta_H_n),
    }
