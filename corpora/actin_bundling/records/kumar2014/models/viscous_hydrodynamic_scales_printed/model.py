import math


def run(inputs, settings):
    d = int(round(inputs["d"]))
    eta = float(inputs["eta"])
    eta_v = float(inputs["eta_v"])
    gamma = float(inputs["gamma"])
    ell = math.sqrt(eta / gamma)
    if d == 1:
        nu = 1.0
    elif d == 2:
        nu = 1.0 + eta_v / eta
    else:
        raise ValueError(f"Unsupported spatial dimension d={d}")
    return {"ell": ell, "nu": nu}
