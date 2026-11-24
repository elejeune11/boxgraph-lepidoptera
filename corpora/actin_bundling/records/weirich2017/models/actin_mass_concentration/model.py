def run(inputs, settings):
    # uM = 1e-6 mol/L; (mol/L) * (g/mol) = g/L = mg/mL
    c = inputs["c_A"] * 1e-6 * inputs["M_w"]
    return {"c": c}


def symbolic():
    import sympy as sp

    c_A, M_w = sp.symbols("c_A M_w", positive=True)
    return {"c": c_A * sp.Rational(1, 10**6) * M_w}
