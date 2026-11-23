import math


def run(inputs, settings):
    xs = inputs["profile"]["x"]
    ys = inputs["profile"]["y"]
    L = inputs["L"]
    num = 0.0
    den = 0.0
    for x, y in zip(xs, ys):
        c = 0.5 * (1.0 + math.cos(2.0 * math.pi * x / L))
        num += y * c
        den += c * c
    return {"a": float(num / den)}
