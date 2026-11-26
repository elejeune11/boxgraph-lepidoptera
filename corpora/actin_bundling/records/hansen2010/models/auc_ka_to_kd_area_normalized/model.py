def run(inputs, settings):
    # When physical concentration c_V(X) (mean C_vasp over span delta_X) is normalized by its
    # area V_TOT_phys = C_vasp * delta_X so that integral(v_norm dX) = 1, mass action
    # v_1 = K_a_norm * v_0 * a_0 implies K_a_phys = K_a_norm / (C_vasp * delta_X), hence
    # K_d_1 = 1 / K_a_phys = C_vasp * delta_X / K_a_norm.
    return {"K_d_1": inputs["C_vasp"] * inputs["delta_X"] / inputs["K_a"]}


def symbolic():
    import sympy as sp

    C_vasp, delta_X, K_a = sp.symbols("C_vasp delta_X K_a", positive=True)
    return {"K_d_1": C_vasp * delta_X / K_a}
