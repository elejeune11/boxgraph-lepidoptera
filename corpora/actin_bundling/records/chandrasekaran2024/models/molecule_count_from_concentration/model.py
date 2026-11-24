import math

AVOGADRO = 6.02214076e23


def run(inputs, settings):
    # C in nM -> C * 1e-9 mol/L; R_drop in um -> (4/3)*pi*R_drop^3 * 1e-15 L
    vol_L = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3) * 1e-15
    return {"N": inputs["C"] * 1e-9 * AVOGADRO * vol_L}


def symbolic():
    import sympy as sp

    C, R_drop = sp.symbols("C R_drop", positive=True)
    return {
        "N": C
        * sp.Float("1e-24")
        * sp.Float("6.02214076e23")
        * sp.Rational(4, 3)
        * sp.pi
        * R_drop**3
    }
