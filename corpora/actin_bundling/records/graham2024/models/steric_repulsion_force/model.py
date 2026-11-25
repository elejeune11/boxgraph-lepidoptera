def run(inputs, settings):
    d = inputs["d"]
    r = inputs["r"]
    # As printed in Fig. S4d: F_steric = k_steric * (d - r) when overlapping (d < r)
    return {"F_steric": inputs["k_steric"] * (d - r) if d < r else 0.0}


def symbolic():
    import sympy as sp

    k_steric, d, r = sp.symbols("k_steric d r", positive=True)
    return {"F_steric": k_steric * (d - r)}
