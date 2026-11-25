def run(inputs, settings):
    A = float(inputs["A"])
    phi_1 = float(inputs["phi_1"])
    Delta_phi = float(inputs["Delta_phi"])
    Delta_psi = float(inputs["Delta_psi"])
    eps = float(inputs["epsilon"])
    f_pp = A * (3.0 * phi_1**2 - 1.0)
    return {"beta_over_Rs": (f_pp / (Delta_phi * Delta_psi)) * eps}


def symbolic():
    import sympy as sp

    A, phi_1, Delta_phi, Delta_psi, epsilon = sp.symbols(
        "A phi_1 Delta_phi Delta_psi epsilon", real=True
    )
    f_pp = A * (3 * phi_1**2 - 1)
    return {"beta_over_Rs": (f_pp / (Delta_phi * Delta_psi)) * epsilon}
