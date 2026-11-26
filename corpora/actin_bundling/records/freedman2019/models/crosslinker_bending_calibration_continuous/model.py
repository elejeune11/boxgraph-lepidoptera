def run(inputs, settings):
    l_long = inputs["l_long"]
    l_short = inputs["l_short"]
    ramp = 22.0 * (l_long / l_short) - 21.0
    ratio = min(4.0, ramp)
    return {
        "stiffness_ratio": ratio,
        "kappa_long_B": ratio * inputs["kappa_short_B"],
    }
