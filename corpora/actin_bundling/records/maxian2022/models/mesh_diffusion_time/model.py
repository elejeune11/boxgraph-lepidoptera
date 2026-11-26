def run(inputs, settings):
    return {"tau_m": (inputs["ell_m"] ** 2) / inputs["six_D_t"]}


def symbolic():
    import sympy as sp

    ell_m, six_D_t = sp.symbols("ell_m six_D_t", positive=True)
    return {"tau_m": ell_m**2 / six_D_t}
