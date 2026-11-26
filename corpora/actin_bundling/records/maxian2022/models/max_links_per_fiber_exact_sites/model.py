import math


def run(inputs, settings):
    cap_per_site = math.ceil(inputs["dsu"] / inputs["cw"])
    n_sites = inputs["L"] / inputs["dsu"] + 1.0
    return {"n_max": cap_per_site * n_sites}
