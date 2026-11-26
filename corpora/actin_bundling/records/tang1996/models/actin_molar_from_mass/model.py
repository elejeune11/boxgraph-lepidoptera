def run(inputs, settings):
    # mg/ml == g/L; (g/L) / (g/mol) = M; * 1e6 -> uM
    return {"c_molar": (inputs["c_mass"] / inputs["M_actin"]) * 1.0e6}


def symbolic():
    import sympy as sp

    c_mass, M_actin = sp.symbols("c_mass M_actin", positive=True)
    return {"c_molar": sp.Integer(1000000) * c_mass / M_actin}
