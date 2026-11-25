def run(inputs, settings):
    rho_A = inputs["rho_A"]
    ell = inputs["ell"]
    ell_A = inputs["ell_A"]
    dt = inputs["dt"]
    stride = inputs["stride"]
    r_A = inputs["r_A"]
    tau_A = 1.0 / r_A
    return {
        "N_A": rho_A * (ell ** 2) / ell_A,
        "rho_A_2x": 2.0 * rho_A,
        "rho_A_3x": 3.0 * rho_A,
        "f_out": 1.0 / (stride * dt),
        "tau_A": tau_A,
        "t_norm": inputs["t_eval"] / tau_A,
        "T": inputs["kT"] / inputs["k_B"],
    }

def symbolic():
    import sympy as sp
    rho_A, ell, ell_A, dt, stride, r_A, t_eval, kT, k_B = sp.symbols(
        "rho_A ell ell_A dt stride r_A t_eval kT k_B", positive=True
    )
    tau_A = 1 / r_A
    return {
        "N_A": rho_A * ell**2 / ell_A,
        "rho_A_2x": 2 * rho_A,
        "rho_A_3x": 3 * rho_A,
        "f_out": 1 / (stride * dt),
        "tau_A": tau_A,
        "t_norm": t_eval / tau_A,
        "T": kT / k_B,
    }
