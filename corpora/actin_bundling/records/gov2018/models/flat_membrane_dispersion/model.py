import math


def run(inputs, settings):
    q = inputs["q"]
    kappa = inputs["kappa"]
    sigma = inputs["sigma"]
    kBT = inputs["kBT"]
    Lambda = inputs["Lambda"]
    H0 = inputs["H0"]
    ns = inputs["ns"]
    n0 = inputs["n0"]
    J = inputs["J"]
    xi = inputs["xi"]
    A = inputs["A"]

    q2 = q * q
    q4 = q2 * q2
    D = Lambda * kBT

    Gh = (kappa * q4 + sigma * q2) / xi
    Gn = (
        D * q2
        + (kappa * Lambda * (H0 ** 2) * n0 / ns) * q2
        - (2.0 * J * Lambda * n0 / ns) * q2 * (1.0 - q2 / ns)
    )
    m12 = (A - kappa * H0 * q2) / xi
    m21 = -(kappa * Lambda * H0 * n0 / ns) * q4

    tr = -(Gh + Gn)
    det = Gh * Gn - m12 * m21
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        omega_real = 0.5 * (tr + math.sqrt(disc))
        omega_imag = 0.0
    else:
        omega_real = 0.5 * tr
        omega_imag = 0.5 * math.sqrt(-disc)
    return {"omega_real": omega_real, "omega_imag": omega_imag}


def symbolic():
    import sympy as sp

    q, kappa, sigma, kBT, Lambda, H0, ns, n0, J, xi, A = sp.symbols(
        "q kappa sigma kBT Lambda H0 ns n0 J xi A", real=True
    )
    q2 = q ** 2
    q4 = q ** 4
    D = Lambda * kBT
    Gh = (kappa * q4 + sigma * q2) / xi
    Gn = (
        D * q2
        + (kappa * Lambda * H0 ** 2 * n0 / ns) * q2
        - (2 * J * Lambda * n0 / ns) * q2 * (1 - q2 / ns)
    )
    m12 = (A - kappa * H0 * q2) / xi
    m21 = -(kappa * Lambda * H0 * n0 / ns) * q4
    tr = -(Gh + Gn)
    det = Gh * Gn - m12 * m21
    return {"omega_real": (tr + sp.sqrt(tr ** 2 - 4 * det)) / 2}
