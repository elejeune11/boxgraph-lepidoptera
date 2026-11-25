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

    if alpha * alpha < chi * chi:
        return {"phi2_bar": float("nan")}
    s = math.sqrt(alpha * alpha - chi * chi) / 6.0
    dx = p1 - p1m
    if branch in ("top", "bot"):
        # Vertical conjugate hyperbola: (p2 - p2m)^2 - (p1 - p1m)^2 = (c2^2 - c1^2)/3 + s
        dy2 = dx * dx + (c2 * c2 - c1 * c1) / 3.0 + s
        if dy2 < 0.0:
            return {"phi2_bar": float("nan")}
        dy = math.sqrt(dy2)
        return {"phi2_bar": p2m + dy if branch == "top" else p2m - dy}
    # Horizontal conjugate hyperbola: (p1 - p1m)^2 - (p2 - p2m)^2 = (c1^2 - c2^2)/3 + s
    off_p = (c1 * c1 - c2 * c2) / 3.0 + s
    dy2 = dx * dx - off_p
    if dy2 < 0.0:
        return {"phi2_bar": float("nan")}
    dy = math.sqrt(dy2)
    return {"phi2_bar": p2m + dy if branch == "right_upper" else p2m - dy}
