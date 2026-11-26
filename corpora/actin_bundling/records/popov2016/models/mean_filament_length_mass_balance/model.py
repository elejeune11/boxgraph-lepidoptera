# 1 uM * 1 um^3 = 602.214076 molecules
MOL_PER_UM_UM3 = 602.214076

def run(inputs, settings):
    n_f_actin = (inputs["A_tot"] - inputs["A_ss"]) * inputs["V"] * MOL_PER_UM_UM3
    return {"L_mean": (n_f_actin / inputs["N_fil"]) * inputs["l_mon"]}

def symbolic():
    import sympy as sp
    A_tot, A_ss, V, N_fil, l_mon = sp.symbols("A_tot A_ss V N_fil l_mon", positive=True)
    return {"L_mean": ((A_tot - A_ss) * V * MOL_PER_UM_UM3 / N_fil) * l_mon}
