def run(inputs, settings):
    # 1 mg/ml = 1 g/L; dividing by M_actin [g/mol] gives mol/L (M); multiply by 1e6 for uM
    return {"c_a_molar": (inputs["c_a"] / inputs["M_actin"]) * 1e6}


def symbolic():
    import sympy as sp

    c_a, M_actin = sp.symbols("c_a M_actin", positive=True)
    return {"c_a_molar": (c_a / M_actin) * sp.Integer(10) ** 6}
