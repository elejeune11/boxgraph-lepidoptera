"""Axial z-scan step count N_steps = depth / dz (p. 6, SI p. 10).

depth is in um (1 um = 1000 nm) and dz is in nm.
"""


def run(inputs, settings):
    return {"N_steps": float(inputs["depth"]) * 1000.0 / float(inputs["dz"])}


def symbolic():
    import sympy as sp

    depth, dz = sp.symbols("depth dz", positive=True)
    return {"N_steps": 1000 * depth / dz}
