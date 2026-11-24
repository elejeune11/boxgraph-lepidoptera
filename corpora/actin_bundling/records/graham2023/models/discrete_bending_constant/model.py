def run(inputs, settings):
    # L_p in um -> 1000 nm; L_cyl in nm; kBT in pN*nm -> k_bend in pN*nm
    return {"k_bend": inputs["N_fil"] * (inputs["L_p"] * 1000.0) * inputs["kBT"] / inputs["L_cyl"]}


def symbolic():
    import sympy as sp

    N_fil, L_p, kBT, L_cyl = sp.symbols("N_fil L_p kBT L_cyl", positive=True)
    return {"k_bend": N_fil * 1000 * L_p * kBT / L_cyl}
