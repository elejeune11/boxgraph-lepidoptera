def run(inputs, settings):
    # Printed verbatim without factor 1/2 (unlike Eqs. 2 and 9).
    dd = inputs["d_ij"] - inputs["d_0"]
    return {"U_branch_str": inputs["K_branch_str"] * dd * dd}

def symbolic():
    import sympy as sp
    K_branch_str, d_ij, d_0 = sp.symbols("K_branch_str d_ij d_0", real=True)
    return {"U_branch_str": K_branch_str * (d_ij - d_0) ** 2}
