def run(inputs, settings):
    A = float(inputs["A"])
    phi_1 = float(inputs["phi_1"])
    Delta_phi = float(inputs["Delta_phi"])
    Delta_psi = float(inputs["Delta_psi"])
    Delta_g = float(inputs["Delta_g"])
    sigma = float(inputs["sigma"])
    d = float(inputs["d"])
    eps = float(inputs["epsilon"])
    R = float(inputs["R"])

    phi_s = phi_1 + eps
    fprime_s = A * (phi_s**3 - phi_s)
    beta = (d - 1.0) * sigma / Delta_psi
    denom = fprime_s * Delta_psi - Delta_g
    R_s = Delta_psi * beta / denom
    beta_over_Rs = denom / Delta_psi
    R_dot = (beta / (R * Delta_phi)) * (1.0 / R_s - 1.0 / R)
    return {
        "beta": beta,
        "R_s": R_s,
        "beta_over_Rs": beta_over_Rs,
        "R_dot": R_dot,
    }


def symbolic():
    import sympy as sp

    A, phi_1, Delta_phi, Delta_psi, Delta_g, sigma, d, epsilon, R = sp.symbols(
        "A phi_1 Delta_phi Delta_psi Delta_g sigma d epsilon R", real=True
    )
    phi_s = phi_1 + epsilon
    fprime_s = A * (phi_s**3 - phi_s)
    beta = (d - 1) * sigma / Delta_psi
    denom = fprime_s * Delta_psi - Delta_g
    R_s = Delta_psi * beta / denom
    beta_over_Rs = denom / Delta_psi
    R_dot = (beta / (R * Delta_phi)) * (1 / R_s - 1 / R)
    return {
        "beta": beta,
        "R_s": R_s,
        "beta_over_Rs": beta_over_Rs,
        "R_dot": R_dot,
    }
