def run(inputs, settings):
    return {
        "k_walk_0": inputs["k_NMII_bind"]
        * (inputs["N_total"] - inputs["N_bound_0"])
        / inputs["N_bound_0"]
    }


def symbolic():
    import sympy as sp
    k_NMII_bind, N_total, N_bound_0 = sp.symbols("k_NMII_bind N_total N_bound_0", positive=True)
    return {"k_walk_0": k_NMII_bind * (N_total - N_bound_0) / N_bound_0}
