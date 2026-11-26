import math


def run(inputs, settings):
    q10 = float(inputs["q_10"])  # in 1/angstrom
    # D in angstrom = 4*pi / (sqrt(3) * q10); convert angstrom -> nm by dividing by 10
    D_nm = (4.0 * math.pi / (math.sqrt(3.0) * q10)) / 10.0
    q11 = math.sqrt(3.0) * q10
    return {"D": D_nm, "q_11": q11}


def symbolic():
    import sympy as sp

    q_10 = sp.symbols("q_10", positive=True)
    return {
        "D": 4 * sp.pi / (10 * sp.sqrt(3) * q_10),
        "q_11": sp.sqrt(3) * q_10,
    }
