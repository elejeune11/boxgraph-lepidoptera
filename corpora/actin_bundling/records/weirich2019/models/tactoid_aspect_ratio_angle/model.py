import math


def run(inputs, settings):
    R = inputs["R"]
    alpha = inputs["alpha"]
    L_D = R * math.sin(alpha)
    L_y = R * (1.0 - math.cos(alpha))
    AR = 1.0 / math.tan(0.5 * alpha)
    return {"AR": AR, "L_D": L_D, "L_y": L_y}


def symbolic():
    import sympy as sp
    R, alpha = sp.symbols("R alpha", positive=True)
    return {
        "AR": sp.cot(alpha / 2),
        "L_D": R * sp.sin(alpha),
        "L_y": R * (1 - sp.cos(alpha)),
    }
