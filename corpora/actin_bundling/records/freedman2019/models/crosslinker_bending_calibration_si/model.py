def run(inputs, settings):
    l_long = inputs["l_long"]
    l_short = inputs["l_short"]
    if l_long <= 200.0:
        ratio = 22.0 * (l_long / l_short) - 21.0
    elif l_long >= 300.0:
        ratio = 4.0
    else:
        ratio = float("nan")
    return {
        "stiffness_ratio": ratio,
        "kappa_long_B": ratio * inputs["kappa_short_B"],
    }
