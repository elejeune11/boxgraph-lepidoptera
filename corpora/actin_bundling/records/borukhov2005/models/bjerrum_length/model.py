# Bjerrum length l_B = e^2 / (eps * k_B * T) in Gaussian CGS units, output in Angstrom.
# CODATA 2018 exact/derived constants in CGS:
#   e = 4.803204712570263e-10 statC
#   k_B = 1.380649e-16 erg/K
#   1 cm = 1e8 Angstrom
_E_CGS = 4.803204712570263e-10
_KB_CGS = 1.380649e-16
_CM_TO_ANGSTROM = 1.0e8
_PREF = (_E_CGS ** 2 / _KB_CGS) * _CM_TO_ANGSTROM


def run(inputs, settings):
    l_B = _PREF / (inputs["eps"] * inputs["T"])
    return {"l_B": float(l_B)}


def symbolic():
    import sympy as sp

    eps, T = sp.symbols("eps T", positive=True)
    return {"l_B": _PREF / (eps * T)}
