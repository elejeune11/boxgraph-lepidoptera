import math


def run(inputs, settings):
    # Omega = 2*pi / twist_pitch (1/nm); l_p in um -> 1000 * l_p in nm; l ~ d
    omega = 2.0 * math.pi / inputs["twist_pitch"]
    P = inputs["dE_kBT"] / (inputs["d"] * (1000.0 * inputs["l_p"]) * (omega**3))
    return {"P": P}


def symbolic():
    import sympy as sp
    dE_kBT, d, l_p, twist_pitch = sp.symbols("dE_kBT d l_p twist_pitch", positive=True)
    omega = 2 * sp.pi / twist_pitch
    return {"P": dE_kBT / (d * 1000 * l_p * omega**3)}
