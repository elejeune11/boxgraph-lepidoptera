import math


def run(inputs, settings):
    cap_per_site = math.ceil(inputs["dsu"] / inputs["cw"])
    return {"n_max": cap_per_site * (inputs["L"] / inputs["dsu"])}
