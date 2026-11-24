import math


def run(inputs, settings):
    M = float(inputs["M"])
    Gamma = float(inputs["Gamma"])
    lam_t = float(inputs["lambda_tilde"])
    alpha = float(inputs["alpha"])
    L1 = float(inputs["L1"])
    L3 = float(inputs["L3"])

    l_c = math.sqrt(lam_t)
    lc2 = lam_t
    t_Q = 1.0 / Gamma
    t_phi = lam_t / M
    gamma = t_Q / t_phi
    E1 = L1 / lc2
    E3 = L3 / lc2
    chi = alpha / lc2
    # Dimensional Eq. (S17) has -(lambda_tilde + alpha/2) grad_bar^2 phi = -(1 + chi/2) grad^2 phi
    lam = 1.0 + 0.5 * chi
    return {
        "l_c": l_c,
        "t_Q": t_Q,
        "t_phi": t_phi,
        "gamma": gamma,
        "E1": E1,
        "E3": E3,
        "chi": chi,
        "lam": lam,
    }


def symbolic():
    import sympy as sp

    M, Gamma, lambda_tilde, alpha, L1, L3 = sp.symbols(
        "M Gamma lambda_tilde alpha L1 L3", positive=True
    )
    l_c = sp.sqrt(lambda_tilde)
    t_Q = 1 / Gamma
    t_phi = lambda_tilde / M
    chi = alpha / lambda_tilde
    return {
        "l_c": l_c,
        "t_Q": t_Q,
        "t_phi": t_phi,
        "gamma": M / (lambda_tilde * Gamma),
        "E1": L1 / lambda_tilde,
        "E3": L3 / lambda_tilde,
        "chi": chi,
        "lam": 1 + chi / 2,
    }
