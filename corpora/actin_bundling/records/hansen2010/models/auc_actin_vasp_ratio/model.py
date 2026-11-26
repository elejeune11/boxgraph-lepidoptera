def run(inputs, settings):
    return {"R_av": inputs["C_actin"] / inputs["C_vasp"]}


def symbolic():
    import sympy as sp

    C_actin, C_vasp = sp.symbols("C_actin C_vasp", positive=True)
    return {"R_av": C_actin / C_vasp}
