import math


def run(inputs, settings):
    p1 = float(inputs["phi1_bar"])
    p1m = float(inputs["phi1_mid"])
    p2m = float(inputs["phi2_mid"])
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    branch = str(settings.get("branch", "upper"))

    r2 = (c1 * c1 + c2 * c2) / 3.0
    dx = p1 - p1m
    dy2 = r2 - dx * dx
    if dy2 < 0.0:
        return {"phi2_bar": float("nan")}
    dy = math.sqrt(dy2)
    return {"phi2_bar": p2m + dy if branch == "upper" else p2m - dy}
