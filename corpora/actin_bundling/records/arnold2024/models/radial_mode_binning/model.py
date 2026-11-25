import math


def run(inputs, settings):
    # Nearest integer rounding (half up, matching standard MATLAB round)
    n_bin = float(math.floor(inputs["n_raw"] + 0.5))
    return {"n_bin": n_bin}
