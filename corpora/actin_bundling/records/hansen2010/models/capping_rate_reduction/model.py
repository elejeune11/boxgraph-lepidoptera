def run(inputs, settings):
    return {"fold_reduction": inputs["K_plus_no_vasp"] / inputs["K_plus_vasp"]}


def symbolic():
    import sympy as sp

    K_plus_no_vasp, K_plus_vasp = sp.symbols("K_plus_no_vasp K_plus_vasp", positive=True)
    return {"fold_reduction": K_plus_no_vasp / K_plus_vasp}
