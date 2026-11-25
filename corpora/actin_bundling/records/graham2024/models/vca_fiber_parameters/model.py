def run(inputs, settings):
    k_bind_fiber = inputs["k_bind"] * inputs["N_VCA"]
    L_Gmax = (inputs["l_mono"] * 1e-3) * inputs["factor_Gmax"] * inputs["N_VCA"]
    return {"k_bind_fiber": k_bind_fiber, "L_Gmax": L_Gmax}


def symbolic():
    import sympy as sp

    k_bind, N_VCA, l_mono, factor_Gmax = sp.symbols(
        "k_bind N_VCA l_mono factor_Gmax", positive=True
    )
    return {
        "k_bind_fiber": k_bind * N_VCA,
        "L_Gmax": l_mono * sp.Rational(1, 1000) * factor_Gmax * N_VCA,
    }
