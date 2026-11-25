import math
from scipy.integrate import solve_ivp


def run(inputs, settings):
    R_f = inputs["R_f"]
    R_b = inputs["R_b"]
    theta_rad = math.radians(inputs["theta"])
    dp_s = inputs["DeltaP_over_sigma"]
    max_step = float(settings.get("max_step", 0.1))
    rtol = float(settings.get("rtol", 1.0e-8))
    atol = float(settings.get("atol", 1.0e-10))

    # At contact line y = R_f, normal vector angle with fiber axis is gamma = pi/2 - theta
    gamma0 = 0.5 * math.pi - theta_rad

    def rhs(y, g):
        return [(dp_s - math.sin(g[0]) / y) / math.cos(g[0])]

    sol = solve_ivp(rhs, [R_f, R_b], [gamma0], method="RK45", max_step=max_step, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"ODE integration failed: {sol.message}")
    # At cut-off y = R_b, gamma(R_b) equals theta_b (angle between contour tangent and y-axis)
    return {"theta_b": math.degrees(sol.y[0, -1])}
