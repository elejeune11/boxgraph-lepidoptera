R_CGS = 8.314462618e7  # erg / (mol * K) = g * cm^2 / (s^2 * mol * K)


def run(inputs, settings):
    M_g_mol = inputs["M"] * 1000.0  # 1 kDa = 1000 g/mol
    omega = inputs["omega"]  # rad/s
    vbar = inputs["vbar"]  # cm^3/g
    rho = inputs["rho"]  # g/cm^3
    T = inputs["T"]  # K
    sigma_o = M_g_mol * (omega ** 2) * (1.0 - vbar * rho) / (R_CGS * T)
    return {"sigma_o": sigma_o}


def symbolic():
    import sympy as sp

    M, omega, vbar, rho, T = sp.symbols("M omega vbar rho T", positive=True)
    return {"sigma_o": (M * 1000) * omega**2 * (1 - vbar * rho) / (R_CGS * T)}
