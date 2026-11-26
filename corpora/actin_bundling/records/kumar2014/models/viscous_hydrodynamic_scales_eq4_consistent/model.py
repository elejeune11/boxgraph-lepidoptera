import math


def run(inputs, settings):
    d = int(round(inputs["d"]))
    eta = float(inputs["eta"])
    eta_v = float(inputs["eta_v"])
    gamma = float(inputs["gamma"])
    if d == 1:
        ell = math.sqrt(eta_v / gamma)
        nu = 1.0
    elif d == 2:
        ell = math.sqrt(eta / gamma)
        nu = 1.0 + eta_v / eta
    else:
        raise ValueError(f"Unsupported spatial dimension d={d}")
    return {"ell": ell, "nu": nu}
