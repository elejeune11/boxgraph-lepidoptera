def run(inputs, settings):
    return {"C_actin": inputs["C_vasp_mon"] / inputs["r_vasp_to_actin"]}


def symbolic():
    import sympy as sp

    C_vasp_mon, r_vasp_to_actin = sp.symbols(
        "C_vasp_mon r_vasp_to_actin", positive=True
    )
    return {"C_actin": C_vasp_mon / r_vasp_to_actin}
