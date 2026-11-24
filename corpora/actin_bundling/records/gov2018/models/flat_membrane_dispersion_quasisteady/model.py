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
    A = inputs["A"]

    q2 = q * q
    omega = -(Lambda * q2 / (ns * (kappa * q2 + sigma))) * (
        A * H0 * kappa * n0
        + sigma * (kBT * ns + kappa * (H0 ** 2) * n0)
        + kappa * kBT * ns * q2
    ) + (2.0 * J * Lambda * n0 * q2 / ns) * (1.0 - q2 / ns)
    return {"omega_real": omega}


def symbolic():
    import sympy as sp

    q, kappa, sigma, kBT, Lambda, H0, ns, n0, J, A = sp.symbols(
        "q kappa sigma kBT Lambda H0 ns n0 J A", real=True
    )
    q2 = q ** 2
    omega = -(Lambda * q2 / (ns * (kappa * q2 + sigma))) * (
        A * H0 * kappa * n0
        + sigma * (kBT * ns + kappa * H0 ** 2 * n0)
        + kappa * kBT * ns * q2
    ) + (2 * J * Lambda * n0 * q2 / ns) * (1 - q2 / ns)
    return {"omega_real": omega}
