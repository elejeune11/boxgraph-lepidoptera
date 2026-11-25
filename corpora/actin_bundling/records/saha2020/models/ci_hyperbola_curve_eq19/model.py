import math


def run(inputs, settings):
    p1 = float(inputs["phi1_bar"])
    p1m = float(inputs["phi1_mid"])
    p2m = float(inputs["phi2_mid"])
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    chi = float(inputs["chi"])
    alpha = float(inputs["alpha"])
    branch = str(settings.get("branch", "right_upper"))

    # Printed Eq. (19): (p1 - p1m)^2 - (p2 - p2m)^2 = offset
    offset = (c1 * c1 + c2 * c2 - chi * chi + alpha * alpha) / 3.0
    dx = p1 - p1m
    dy2 = dx * dx - offset
    if dy2 < 0.0:
        return {"phi2_bar": float("nan")}
    dy = math.sqrt(dy2)
    if branch == "bot":
        return {"phi2_bar": p2m - dy}
    return {"phi2_bar": p2m + dy}
