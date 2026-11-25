from scipy.optimize import brentq


def run(inputs, settings):
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Om = float(inputs["Omega"])
    Omd0 = float(inputs["Omega_d0"])
    Omd = float(inputs["Omega_d"])
    Z = float(inputs["Z"])
    B = float(inputs["B"])
    Nbar = float(inputs["Nbar"])

    ratio = A / Kd

    def f(na):
        return Om0 * (1.0 + Om * na * na) * (Nbar - na) - (Omd0 + Omd * ratio * na) * na

    na_h = float(brentq(f, 0.0, Nbar, xtol=1e-14, rtol=1e-14))
    ni_h = Nbar - na_h
    c_h = ratio * na_h
    nh_minus = 2.0 * na_h - Nbar
    pi_prime_h = 2.0 * Z * c_h - 3.0 * B * c_h * c_h
    return {
        "C_H": c_h,
        "Na_H": na_h,
        "Ni_H": ni_h,
        "NH_minus": nh_minus,
        "Pi_prime_H": pi_prime_h,
    }
