def run(inputs, settings):
    th0 = -360.0 * inputs["s_0"] / inputs["p_0"]
    th1 = -360.0 * inputs["s_1"] / inputs["p_1"]
    return {
        "theta_0": th0,
        "theta_1": th1,
        "Delta_theta": th0 - th1,
    }


def symbolic():
    import sympy as sp

    p_0, s_0, p_1, s_1 = sp.symbols("p_0 s_0 p_1 s_1", positive=True)
    return {
        "theta_0": -360 * s_0 / p_0,
        "theta_1": -360 * s_1 / p_1,
        "Delta_theta": 360 * (s_1 / p_1 - s_0 / p_0),
    }
