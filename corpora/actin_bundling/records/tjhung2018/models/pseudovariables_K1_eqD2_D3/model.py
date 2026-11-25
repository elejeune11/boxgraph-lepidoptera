import math


def run(inputs, settings):
    K = float(inputs["K"])
    K_1 = float(inputs["K_1"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    phi = float(inputs["phi"])

    nu = lam - K_1 - 0.5 * zeta
    if abs(K_1) < 1e-10:
        # Smooth K_1 -> 0 limit: psi_bar -> exp(-2*nu*phi/K) = exp((zeta - 2*lam)*phi/K)
        expo = -2.0 * nu * phi / K
        psi_bar = math.exp(expo)
    else:
        base = 1.0 + 2.0 * K_1 * phi / K
        if base <= 0.0:
            raise ValueError(f"1 + 2*K_1*phi/K = {base} <= 0; K(phi) must remain positive")
        psi_bar = base ** (-nu / K_1)

    if abs(nu) < 1e-10:
        # Smooth nu -> 0 limit: psi_{K_1} -> (K/(2*K_1))*ln(1 + 2*K_1*phi/K) or phi as K_1 -> 0
        if abs(K_1) < 1e-10:
            psi_K1 = phi
        else:
            psi_K1 = (K / (2.0 * K_1)) * math.log(1.0 + 2.0 * K_1 * phi / K)
    else:
        psi_K1 = -(K / (2.0 * nu)) * (psi_bar - 1.0)

    return {"nu": nu, "psi_bar": psi_bar, "psi_K1": psi_K1}


def symbolic():
    import sympy as sp

    K, K_1, zeta, lam, phi = sp.symbols("K K_1 zeta lam phi", real=True)
    nu = lam - K_1 - zeta / 2
    psi_bar = (1 + 2 * K_1 * phi / K) ** (-nu / K_1)
    psi_K1 = -(K / (2 * nu)) * (psi_bar - 1)
    return {"nu": nu, "psi_bar": psi_bar, "psi_K1": psi_K1}
