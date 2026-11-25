import math


def _psi_and_g(A, K, alpha, phi):
    if abs(alpha) < 1e-3:
        # Smooth Taylor series avoiding 1/alpha^4 catastrophic cancellation near alpha = 0
        # psi(phi) = sum_{n=0}^inf alpha^n phi^(n+1) / (n+1)!
        # g(phi) = A * sum_{n=0}^inf (alpha^n / n!) * (phi^(n+4)/(n+4) - phi^(n+2)/(n+2))
        psi_val = 0.0
        g_val = 0.0
        fact = 1.0
        apow = 1.0
        for n in range(14):
            if n > 0:
                fact *= n
                apow *= alpha
            psi_val += apow * (phi ** (n + 1)) / (fact * (n + 1))
            g_val += A * (apow / fact) * ((phi ** (n + 4)) / (n + 4) - (phi ** (n + 2)) / (n + 2))
        return psi_val, g_val
    e = math.exp(alpha * phi)
    psi_val = (e - 1.0) / alpha
    bracket = (
        -(alpha**3) * phi * (phi**2 - 1.0)
        + 6.0
        - 6.0 * alpha * phi
        + (alpha**2) * (3.0 * phi**2 - 1.0)
    )
    g_val = A * (6.0 / alpha**4 - 1.0 / alpha**2 - (e / alpha**4) * bracket)
    return psi_val, g_val


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    phi = float(inputs["phi"])
    alpha = (zeta - 2.0 * lam) / K
    psi_val, g_val = _psi_and_g(A, K, alpha, phi)
    fprime = A * (phi**3 - phi)
    P_val = psi_val * fprime - g_val
    return {"psi": psi_val, "g": g_val, "P": P_val}


def symbolic():
    import sympy as sp

    A, K, zeta, lam, phi = sp.symbols("A K zeta lam phi", real=True)
    d = zeta - 2 * lam
    e = sp.exp(d * phi / K)
    psi_val = (K / d) * (e - 1)
    bracket = -(d**3) * phi * (phi**2 - 1) + 6 * K**3 - 6 * d * K**2 * phi + d**2 * K * (3 * phi**2 - 1)
    g_val = 6 * A * K**4 / d**4 - A * K**2 / d**2 - (A * K * e / d**4) * bracket
    P_val = psi_val * A * (phi**3 - phi) - g_val
    return {"psi": psi_val, "g": g_val, "P": P_val}
