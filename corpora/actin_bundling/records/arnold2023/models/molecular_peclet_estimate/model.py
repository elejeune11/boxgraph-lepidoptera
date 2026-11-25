def run(inputs, settings):
    # Convert l from nm to um (1 nm = 1e-3 um)
    l_um = inputs["l"] * 1e-3
    pe_mol = inputs["v"] * l_um / inputs["D_lip"]
    return {"Pe_mol": float(pe_mol)}


def symbolic():
    import sympy as sp
    v, l, D_lip = sp.symbols("v l D_lip", positive=True)
    return {"Pe_mol": v * (l * sp.Rational(1, 1000)) / D_lip}
