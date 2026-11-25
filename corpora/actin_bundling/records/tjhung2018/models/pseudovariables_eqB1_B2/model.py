import math


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    phi = float(inputs["phi"])
    d = zeta - 2.0 * lam
    # Printed Eq. (B1) and Eq. (B2) verbatim (note -6*zeta*K^2*phi in Eq. B2)
    if abs(d) < 1e-12:
        if abs(zeta) < 1e-12:
            psi_val = phi
            g_val = A * (-0.5 * phi**2 + 0.25 * phi**4)
        else:
            # Printed Eq. (B2) genuinely diverges as d -> 0 when zeta != 0 (lambda = zeta/2)
            raise ValueError("Printed Eq. (B2) diverges at zeta - 2*lambda = 0 when zeta != 0 due to -6*zeta*K^2*phi")
    else:
        e = math.exp((d / K) * phi)
        psi_val = (K / d) * (e - 1.0)
        bracket = (
            -(d**3) * phi * (phi**2 - 1.0)
            + 6.0 * K**3
            - 6.0 * zeta * K**2 * phi
            + (d**2) * K * (3.0 * phi**2 - 1.0)
        )
        g_val = 6.0 * A * K**4 / d**4 - A * K**2 / d**2 - (A * K * e / d**4) * bracket
    fprime = A * (phi**3 - phi)
    P_val = psi_val * fprime - g_val
    return {"psi": psi_val, "g": g_val, "P": P_val}


def symbolic():
    import sympy as sp

    A, K, zeta, lam, phi = sp.symbols("A K zeta lam phi", real=True)
    d = zeta - 2 * lam
    e = sp.exp(d * phi / K)
    psi_val = (K / d) * (e - 1)
    bracket = -(d**3) * phi * (phi**2 - 1) + 6 * K**3 - 6 * zeta * K**2 * phi + d**2 * K * (3 * phi**2 - 1)
    g_val = 6 * A * K**4 / d**4 - A * K**2 / d**2 - (A * K * e / d**4) * bracket
    P_val = psi_val * A * (phi**3 - phi) - g_val
    return {"psi": psi_val, "g": g_val, "P": P_val}
