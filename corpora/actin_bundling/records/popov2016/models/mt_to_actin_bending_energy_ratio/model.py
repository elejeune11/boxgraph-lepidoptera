def run(inputs, settings):
    # At equal number of monomers m per cylindrical segment, l_cyl = m * l_mon,
    # so eps_MT / eps_actin = (l_p_MT / l_mon_MT) / (l_p_actin / l_mon_actin).
    # At equal physical segment length l_cyl, eps_MT / eps_actin = l_p_MT / l_p_actin.
    ratio_len = inputs["l_p_MT"] / inputs["l_p_actin"]
    ratio_mon = ratio_len * (inputs["l_mon_actin"] / inputs["l_mon_MT"])
    return {"ratio_same_monomers": ratio_mon, "ratio_same_length": ratio_len}

def symbolic():
    import sympy as sp
    l_p_MT, l_mon_MT, l_p_actin, l_mon_actin = sp.symbols(
        "l_p_MT l_mon_MT l_p_actin l_mon_actin", positive=True
    )
    ratio_len = l_p_MT / l_p_actin
    ratio_mon = ratio_len * l_mon_actin / l_mon_MT
    return {"ratio_same_monomers": ratio_mon, "ratio_same_length": ratio_len}
