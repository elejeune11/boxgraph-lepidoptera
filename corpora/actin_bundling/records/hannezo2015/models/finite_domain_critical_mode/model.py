import math


def run(inputs, settings):
    a = inputs["alpha"]
    phi = inputs["phi"]
    x = a * (phi ** 0.25) / math.pi
    fl = math.floor(x)
    if fl >= 1 and fl > (x * x - fl * fl):
        ic = fl
    else:
        ic = max(1, math.ceil(x))
    k = ic * math.pi
    chi_c0 = 1.0 + (k * k) / (a * a) + phi * (1.0 + (a * a) / (k * k))
    return {"i_c": float(ic), "n_rings": float(ic) / 2.0, "chi_c0": chi_c0}
