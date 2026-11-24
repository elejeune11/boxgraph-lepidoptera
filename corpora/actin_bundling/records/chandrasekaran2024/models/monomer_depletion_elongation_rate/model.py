def run(inputs, settings):
    return {
        "v_elong": inputs["k_grow"]
        * inputs["N_fil"]
        * (1.0 - inputs["F_actin"] / inputs["T_actin"])
    }


def symbolic():
    import sympy as sp

    k_grow, N_fil, F_actin, T_actin = sp.symbols(
        "k_grow N_fil F_actin T_actin", positive=True
    )
    return {"v_elong": k_grow * N_fil * (1 - F_actin / T_actin)}
