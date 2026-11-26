def _at_end(curve):
    xs, ys = curve["x"], curve["y"]
    idx = max(range(len(xs)), key=lambda i: xs[i])
    return float(ys[idx])

def run(inputs, settings):
    rgi = float(inputs["Rg_i"])
    pairs = [
        (0.01, "c_0p01"),
        (0.02, "c_0p02"),
        (0.05, "c_0p05"),
        (0.1,  "c_0p1"),
        (0.2,  "c_0p2"),
        (0.5,  "c_0p5"),
    ]
    xs = [ra for ra, _ in pairs]
    ys = [rgi * _at_end(inputs[k]) for _, k in pairs]
    return {"Rg_final": {"x": xs, "y": ys}}
