def run(inputs, settings):
    p1 = float(inputs["phi1"])
    p2 = float(inputs["phi2"])
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    alpha = float(inputs["alpha"])

    mu1_eq = 2.0 * (p1 - c11) * (p1 - c12) * (2.0 * p1 - c11 - c12) + chi * p2 + 2.0 * chip * p1 * p2 * p2
    # Printed verbatim as Eq. (6) on p. 5: first term prints phi1 instead of phi2.
    mu2_eq = 2.0 * (p1 - c21) * (p1 - c22) * (2.0 * p1 - c21 - c22) + chi * p1 + 2.0 * chip * p2 * p1 * p1
    mu1_neq = mu1_eq + alpha * p2
    mu2_neq = mu2_eq - alpha * p1
    return {
        "mu1_eq": mu1_eq,
        "mu2_eq": mu2_eq,
        "mu1_neq": mu1_neq,
        "mu2_neq": mu2_neq,
    }


def symbolic():
    import sympy as sp

    phi1, phi2, c11, c12, c21, c22, chi, chip, alpha = sp.symbols(
        "phi1 phi2 c11 c12 c21 c22 chi chip alpha", real=True
    )
    mu1_eq = 2 * (phi1 - c11) * (phi1 - c12) * (2 * phi1 - c11 - c12) + chi * phi2 + 2 * chip * phi1 * phi2**2
    mu2_eq = 2 * (phi1 - c21) * (phi1 - c22) * (2 * phi1 - c21 - c22) + chi * phi1 + 2 * chip * phi2 * phi1**2
    return {
        "mu1_eq": mu1_eq,
        "mu2_eq": mu2_eq,
        "mu1_neq": mu1_eq + alpha * phi2,
        "mu2_neq": mu2_eq - alpha * phi1,
    }
