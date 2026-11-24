import math


def run(inputs, settings):
    R = inputs["R"]
    D_B = inputs["D_B"]
    psi = inputs["psi"]
    phi_0_B = inputs["phi_0_B"]
    phi_0_A = inputs["phi_0_A"]
    gb = inputs["gamma_beta"]
    k_s = inputs["k_s"]
    k_b = inputs["k_b"]
    Q = inputs["Q"]

    V = (4.0 * math.pi / 3.0) * R**3
    dV_dt = (
        (4.0 * math.pi * D_B / psi) * (phi_0_B * R - gb)
        + V * (k_s * phi_0_A / psi - k_b)
        + Q * phi_0_A / psi
    )
    dR_dt = dV_dt / (4.0 * math.pi * R**2)
    return {"dV_dt": dV_dt, "dR_dt": dR_dt}


def symbolic():
    import sympy as sp

    R, D_B, psi, phi_0_B, phi_0_A, gamma_beta, k_s, k_b, Q = sp.symbols(
        "R D_B psi phi_0_B phi_0_A gamma_beta k_s k_b Q", positive=True
    )
    V = sp.Rational(4, 3) * sp.pi * R**3
    dV_dt = (
        (4 * sp.pi * D_B / psi) * (phi_0_B * R - gamma_beta)
        + V * (k_s * phi_0_A / psi - k_b)
        + Q * phi_0_A / psi
    )
    dR_dt = sp.simplify(dV_dt / (4 * sp.pi * R**2))
    return {"dV_dt": dV_dt, "dR_dt": dR_dt}
