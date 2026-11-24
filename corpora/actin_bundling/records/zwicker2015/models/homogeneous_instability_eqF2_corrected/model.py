import math


def run(inputs, settings):
    D_B = inputs["D_B"]
    w = inputs["w"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    ksum = k_f + k_b
    rad = 1.0 - 2.0 * w * math.sqrt(ksum / D_B)
    if rad <= 0.0:
        hw = 0.0
    else:
        hw = (ksum / math.sqrt(3.0)) * math.sqrt(rad)
    return {
        "half_width_rate": hw,
        "margin": hw - abs(k_f - k_b),
        "max_rate_sum_norm": 0.25,
    }


def symbolic():
    import sympy as sp

    D_B, w, k_f, k_b = sp.symbols("D_B w k_f k_b", positive=True)
    ksum = k_f + k_b
    hw = (ksum / sp.sqrt(3)) * sp.sqrt(1 - 2 * w * sp.sqrt(ksum / D_B))
    return {
        "half_width_rate": hw,
        "margin": hw - sp.Abs(k_f - k_b),
        "max_rate_sum_norm": sp.Rational(1, 4),
    }
