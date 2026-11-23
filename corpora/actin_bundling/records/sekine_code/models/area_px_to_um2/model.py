"""Convert pixel area count (Nanocluster_motion.m line 40: [6 53] px) to physical area (um^2)."""


def run(inputs, settings):
    return {"A_phys": float(inputs["A_px"]) * (float(inputs["dx"]) ** 2)}


def symbolic():
    import sympy as sp

    A_px, dx = sp.symbols("A_px dx", positive=True)
    return {"A_phys": A_px * dx**2}
