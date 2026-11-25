def run(inputs, settings):
    k = inputs["k"]      # pN/um
    d = inputs["d"]      # um
    L_D = inputs["L_D"]  # um
    f = k * d            # pN
    Delta_F_pN_nm = 0.5 * k * (d ** 2) * 1000.0
    F_0_pN_nm = 0.5 * k * (L_D ** 2) * 1000.0
    return {"f": f, "Delta_F": Delta_F_pN_nm, "F_0": F_0_pN_nm}


def symbolic():
    import sympy as sp
    k, d, L_D = sp.symbols("k d L_D", positive=True)
    return {
        "f": k * d,
        "Delta_F": 500 * k * d ** 2,
        "F_0": 500 * k * L_D ** 2,
    }
