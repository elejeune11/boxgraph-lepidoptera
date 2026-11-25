"""Longitudinal scale length extrapolation to the wing surface plane (SI p. 11, Fig. S9C).

In the longitudinal cross-section perpendicular to the wing plane, a straight scale inclined
at angle theta (deg) with distal tip height z_tip (um) above the wing plane intersects the
wing plane at distance L = z_tip / sin(theta).
"""

import math


def run(inputs, settings):
    z_tip = float(inputs["z_tip"])
    theta_rad = math.radians(float(inputs["theta"]))
    return {"L": z_tip / math.sin(theta_rad)}


def symbolic():
    import sympy as sp

    z_tip, theta = sp.symbols("z_tip theta", positive=True)
    return {"L": z_tip / sp.sin(theta * sp.pi / 180)}
