def run(inputs, settings):
    return {"C_vasp_tet": inputs["C_vasp_mon"] / inputs["valency"]}


def symbolic():
    import sympy as sp

    C_vasp_mon, valency = sp.symbols("C_vasp_mon valency", positive=True)
    return {"C_vasp_tet": C_vasp_mon / valency}
