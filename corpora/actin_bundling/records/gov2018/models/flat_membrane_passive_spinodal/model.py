def run(inputs, settings):
    kBT = inputs["kBT"]
    ns = inputs["ns"]
    n0 = inputs["n0"]
    kappa = inputs["kappa"]
    H0 = inputs["H0"]
    J_c_tensioned = (kBT * ns + kappa * (H0 ** 2) * n0) / (2.0 * n0)
    J_c_tensionless = (kBT * ns) / (2.0 * n0)
    return {
        "J_c_tensioned": J_c_tensioned,
        "J_c_tensionless": J_c_tensionless,
    }


def symbolic():
    import sympy as sp

    kappa, kBT, H0, ns, n0 = sp.symbols("kappa kBT H0 ns n0", positive=True)
    return {
        "J_c_tensioned": (kBT * ns + kappa * H0 ** 2 * n0) / (2 * n0),
        "J_c_tensionless": (kBT * ns) / (2 * n0),
    }
