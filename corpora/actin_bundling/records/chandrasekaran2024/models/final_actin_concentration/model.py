import math

AVOGADRO = 6.02214076e23


def run(inputs, settings):
    # L_fil in um, dl in nm -> monomers per filament = (L_fil * 1000) / dl
    n_monomers = inputs["N_fil"] * (inputs["L_fil"] * 1000.0) / inputs["dl"]
    # R_drop in um -> volume in liters: (4/3)*pi*(R_drop * 1e-6 m)^3 * 1e3 L/m^3 = (4/3)*pi*R_drop^3 * 1e-15 L
    vol_L = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3) * 1e-15
    # Molar -> uM (* 1e6)
    return {"C_actin": (n_monomers / (AVOGADRO * vol_L)) * 1e6}


def symbolic():
    import sympy as sp

    N_fil, L_fil, dl, R_drop = sp.symbols("N_fil L_fil dl R_drop", positive=True)
    # 1000 / (1e-15) * 1e6 = 1e24
    return {
        "C_actin": (N_fil * L_fil * sp.Float("1e24"))
        / (dl * sp.Float("6.02214076e23") * sp.Rational(4, 3) * sp.pi * R_drop**3)
    }
