def run(inputs, settings):
    D_M_um = inputs["D_M"] * 1e-3  # nm -> um
    L_M_um = inputs["L_M"]         # um
    L_a_um = inputs["L_a"] * 1e-3  # nm -> um
    c_a = inputs["c_a"]            # 1/um^3
    kBT = inputs["kBT"]            # pN*nm
    n_excl = D_M_um * L_M_um * L_a_um * c_a
    return {"F_0": n_excl * kBT, "F_0_kBT": n_excl}


def symbolic():
    import sympy as sp
    D_M, L_M, L_a, c_a, kBT = sp.symbols("D_M L_M L_a c_a kBT", positive=True)
    n_excl = (D_M * sp.Rational(1, 1000)) * L_M * (L_a * sp.Rational(1, 1000)) * c_a
    return {"F_0": n_excl * kBT, "F_0_kBT": n_excl}
