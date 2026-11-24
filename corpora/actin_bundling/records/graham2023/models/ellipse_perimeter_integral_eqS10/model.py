import math


def run(inputs, settings):
    a = inputs["a"]
    b = inputs["b"]
    n_quad = max(4, int(round(settings.get("n_quad", 200))))
    dtheta = (0.5 * math.pi) / n_quad
    total = 0.0
    for k in range(n_quad):
        theta = (k + 0.5) * dtheta
        s = math.sin(theta)
        c = math.cos(theta)
        total += math.sqrt((a * s) ** 2 + (b * c) ** 2)
    return {"L": 4.0 * total * dtheta}
