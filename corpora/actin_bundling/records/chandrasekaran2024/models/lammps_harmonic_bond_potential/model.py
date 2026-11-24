def run(inputs, settings):
    # LAMMPS bond_style harmonic convention: E = K * (r - r0)^2 (factor 1/2 absorbed in K)
    return {"U_bond": inputs["K_bond"] * (inputs["r"] - inputs["r0"]) ** 2}


def symbolic():
    import sympy as sp

    K_bond, r, r0 = sp.symbols("K_bond r r0", positive=True)
    return {"U_bond": K_bond * (r - r0) ** 2}
