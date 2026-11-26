def run(inputs, settings):
    L_x = inputs["L_bundle"] + 2.0 * inputs["d_skin"]
    l_vox_um = inputs["l_vox"] * 1e-3
    N_x = L_x / l_vox_um
    # 1 um^3 = 1e-15 L
    V_r_L = L_x * inputs["L_y"] * inputs["L_z"] * 1e-15
    return {"L_x": L_x, "N_x": N_x, "V_r": V_r_L}


def symbolic():
    import sympy as sp
    L_bundle, d_skin, L_y, L_z, l_vox = sp.symbols("L_bundle d_skin L_y L_z l_vox", positive=True)
    L_x = L_bundle + 2 * d_skin
    return {
        "L_x": L_x,
        "N_x": L_x / (l_vox * sp.Float("1e-3")),
        "V_r": L_x * L_y * L_z * sp.Float("1e-15"),
    }
