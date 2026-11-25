"""Scale presumptive rectangular bounding-box volume and surface area (Fig. 5F, p. 7, SI p. 11)."""


def run(inputs, settings):
    L = float(inputs["L"])
    W = float(inputs["W"])
    T = float(inputs["T"])
    return {
        "V": L * W * T,
        "S": 2.0 * (L * W + L * T + W * T),
    }


def symbolic():
    import sympy as sp

    L, W, T = sp.symbols("L W T", positive=True)
    return {
        "V": L * W * T,
        "S": 2 * (L * W + L * T + W * T),
    }
