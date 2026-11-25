def run(inputs, settings):
    r_um = float(inputs["r_solid"]) * 1e-3
    R_um = float(inputs["R_drop"])
    phi_pct = 100.0 * float(inputs["N"]) * ((r_um / R_um) ** 3)
    return {"phi": phi_pct}


def symbolic():
    import sympy as sp

    N, r_solid, R_drop = sp.symbols("N r_solid R_drop", positive=True)
    return {"phi": 100 * N * (r_solid * sp.Rational(1, 1000) / R_drop) ** 3}
