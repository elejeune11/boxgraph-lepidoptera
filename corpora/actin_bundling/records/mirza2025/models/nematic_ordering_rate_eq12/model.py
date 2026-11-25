def run(inputs, settings):
    # c_eff_energy represents (rho * lambda_o - 2 * a) in discrete energy units [pN*um]
    s_dot = inputs["c_eff_energy"] * inputs["S_0"] / inputs["eta_rot"]
    return {
        "S_dot": s_dot,
        "S_dot_norm": s_dot / inputs["r_A"],
    }


def symbolic():
    import sympy as sp
    S_0, c_eff_energy, eta_rot, r_A = sp.symbols("S_0 c_eff_energy eta_rot r_A", real=True)
    s_dot = c_eff_energy * S_0 / eta_rot
    return {
        "S_dot": s_dot,
        "S_dot_norm": s_dot / r_A,
    }
