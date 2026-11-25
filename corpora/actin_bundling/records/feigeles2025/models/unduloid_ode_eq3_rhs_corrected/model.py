import math
from scipy.integrate import solve_ivp


def _slope_corrected(y, a, c):
    u = a * y - c / y
    if u >= 1.0:
        return 0.0
    return math.sqrt(1.0 / (u * u) - 1.0)


def run(inputs, settings):
    R_f = inputs["R_f"]
    R_b = inputs["R_b"]
    a = inputs["DeltaP_over_2sigma"]
    c = inputs["c"]
    max_step = float(settings.get("max_step", 0.1))
    rtol = float(settings.get("rtol", 1.0e-8))
    atol = float(settings.get("atol", 1.0e-10))

    R_max = (1.0 + math.sqrt(1.0 + 4.0 * a * c)) / (2.0 * a)

    slope_rf = _slope_corrected(R_f, a, c)
    slope_rb = _slope_corrected(R_b, a, c)
    theta_contact = math.degrees(math.atan(slope_rf))
    theta_b_out = math.degrees(math.atan2(1.0, slope_rb))

    def rhs(x, y):
        return [_slope_corrected(y[0], a, c)]

    def hit_rb(x, y):
        return y[0] - R_b

    hit_rb.terminal = True
    hit_rb.direction = 1

    def hit_apex(x, y):
        return y[0] - (R_max - 1.0e-7)

    hit_apex.terminal = True
    hit_apex.direction = 1

    x_upper = max(100.0, 20.0 * R_max)
    sol_b = solve_ivp(rhs, [0.0, x_upper], [R_f], method="RK45", events=hit_rb, max_step=max_step, rtol=rtol, atol=atol)
    if not sol_b.success or len(sol_b.t_events[0]) == 0:
        raise RuntimeError("ODE failed to reach cut-off radius R_b")
    x_b = float(sol_b.t_events[0][0])

    sol_apex = solve_ivp(rhs, [0.0, x_upper], [R_f], method="RK45", events=hit_apex, max_step=max_step, rtol=rtol, atol=atol)
    if not sol_apex.success or len(sol_apex.t_events[0]) == 0:
        raise RuntimeError("ODE failed to reach apex radius R_max")
    x_apex = float(sol_apex.t_events[0][0])

    aspect_ratio_eq = R_max / x_apex

    return {
        "theta_contact": theta_contact,
        "theta_b_out": theta_b_out,
        "x_b": x_b,
        "x_apex": x_apex,
        "aspect_ratio_eq": aspect_ratio_eq,
    }
