def run(inputs, settings):
    return {
        "vasp_to_actin": inputs["C_vasp"] / inputs["C_actin"],
        "actin_to_vasp": inputs["C_actin"] / inputs["C_vasp"],
    }


def symbolic():
    import sympy as sp

    C_vasp, C_actin = sp.symbols("C_vasp C_actin", positive=True)
    return {
        "vasp_to_actin": C_vasp / C_actin,
        "actin_to_vasp": C_actin / C_vasp,
    }
