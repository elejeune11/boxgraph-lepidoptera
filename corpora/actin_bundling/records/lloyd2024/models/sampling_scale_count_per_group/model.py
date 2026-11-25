"""Total wing scale count per treatment group: N_scales = n_ind * n_scales_per_ind (SI p. 6)."""


def run(inputs, settings):
    return {"N_scales": float(inputs["n_ind"]) * float(inputs["n_scales_per_ind"])}


def symbolic():
    import sympy as sp

    n_ind, n_scales_per_ind = sp.symbols("n_ind n_scales_per_ind", positive=True)
    return {"N_scales": n_ind * n_scales_per_ind}
