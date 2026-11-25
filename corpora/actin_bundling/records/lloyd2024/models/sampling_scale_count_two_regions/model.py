"""Total wing scale count in a two-region cover/ground SEM sampling design (p. 12, SI pp. 4-5)."""


def run(inputs, settings):
    return {
        "N_scales": float(inputs["n_ind"])
        * float(inputs["n_regions"])
        * (float(inputs["n_cover"]) + float(inputs["n_ground"]))
    }


def symbolic():
    import sympy as sp

    n_ind, n_regions, n_cover, n_ground = sp.symbols("n_ind n_regions n_cover n_ground", positive=True)
    return {"N_scales": n_ind * n_regions * (n_cover + n_ground)}
