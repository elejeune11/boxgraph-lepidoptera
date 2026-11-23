import math


def run(inputs, settings):
    val = inputs["alpha"] * (inputs["phi"] ** 0.25) / (math.sqrt(2.0) * math.pi)
    return {"i_0": float(math.floor(val))}
