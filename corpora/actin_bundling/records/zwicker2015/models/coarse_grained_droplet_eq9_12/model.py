import math


def run(inputs, settings):
    R = inputs["R"]
    psi = inputs["psi"]
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    k = inputs["k"]
    Q = inputs["Q"]
    phi_0_A = inputs["phi_0_A"]
    phi_0_B = inputs["phi_0_B"]
    phibar = inputs["phibar"]
    V_s = inputs["V_s"]
    N = inputs["N"]

    V = (4.0 * math.pi / 3.0) * R**3
    V_tot = N * V
    J = 4.0 * math.pi * D_B * (gb - R * phi_0_B)
    dV_dt = ((k_f * phi_0_A - k_b * psi + k * phi_0_A * psi) * V + Q * phi_0_A - J) / psi
    dR_dt = dV_dt / (4.0 * math.pi * R**2)
    phi_0_A_cons = phibar - ((V_s - V_tot) / V_s) * phi_0_B - (V_tot / V_s) * psi
    dphi0B_dt = (
        k_f * phi_0_A
        - k_b * phi_0_B
        + k * phi_0_A * phi_0_B
        + (N * J) / (V_s - V_tot)
    )
    return {
        "dV_dt": dV_dt,
        "dR_dt": dR_dt,
        "phi_0_A_cons": phi_0_A_cons,
        "dphi0B_dt": dphi0B_dt,
    }


def symbolic():
    import sympy as sp

    R, psi, D_B, gamma_beta, k_f, k_b, k, Q, phi_0_A, phi_0_B, phibar, V_s, N = sp.symbols(
        "R psi D_B gamma_beta k_f k_b k Q phi_0_A phi_0_B phibar V_s N", positive=True
    )
    V = sp.Rational(4, 3) * sp.pi * R**3
    V_tot = N * V
    J = 4 * sp.pi * D_B * (gamma_beta - R * phi_0_B)
    dV_dt = ((k_f * phi_0_A - k_b * psi + k * phi_0_A * psi) * V + Q * phi_0_A - J) / psi
    dR_dt = sp.simplify(dV_dt / (4 * sp.pi * R**2))
    phi_0_A_cons = phibar - ((V_s - V_tot) / V_s) * phi_0_B - (V_tot / V_s) * psi
    dphi0B_dt = k_f * phi_0_A - k_b * phi_0_B + k * phi_0_A * phi_0_B + (N * J) / (V_s - V_tot)
    return {
        "dV_dt": dV_dt,
        "dR_dt": dR_dt,
        "phi_0_A_cons": phi_0_A_cons,
        "dphi0B_dt": dphi0B_dt,
    }
