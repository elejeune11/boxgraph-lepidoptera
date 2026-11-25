def run(inputs, settings):
    return {"C_tet": inputs["C_mono"] / inputs["valency"]}


def symbolic():
    import sympy as sp

    C_mono, valency = sp.symbols("C_mono valency", positive=True)
    return {"C_tet": C_mono / valency}
