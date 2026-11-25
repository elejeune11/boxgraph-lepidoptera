def run(inputs, settings):
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Om0 = float(inputs["Omega0"])
    Omd0 = float(inputs["Omega_d0"])
    Nbar = float(inputs["Nbar"])

    na_h = Om0 * Nbar / (Om0 + Omd0)
    ni_h = Nbar - na_h
    c_h = (A / Kd) * na_h
    nh_minus = 2.0 * na_h - Nbar
    return {
        "C_H": c_h,
        "Na_H": na_h,
        "Ni_H": ni_h,
        "NH_minus": nh_minus,
    }


def symbolic():
    import sympy as sp

    A, Kd, Omega0, Omega_d0, Nbar = sp.symbols("A Kd Omega0 Omega_d0 Nbar", positive=True)
    na_h = Omega0 * Nbar / (Omega0 + Omega_d0)
    ni_h = Nbar - na_h
    c_h = (A / Kd) * na_h
    nh_minus = 2 * na_h - Nbar
    return {
        "C_H": c_h,
        "Na_H": na_h,
        "Ni_H": ni_h,
        "NH_minus": nh_minus,
    }
