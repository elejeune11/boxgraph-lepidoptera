def run(inputs, settings):
    # L_p in um -> 1000 nm; kBT in pN*nm -> k_b_fil in pN*nm^2
    return {"k_b_fil": inputs["L_p"] * 1000.0 * inputs["kBT"]}


def symbolic():
    import sympy as sp

    L_p, kBT = sp.symbols("L_p kBT", positive=True)
    return {"k_b_fil": 1000 * L_p * kBT}
