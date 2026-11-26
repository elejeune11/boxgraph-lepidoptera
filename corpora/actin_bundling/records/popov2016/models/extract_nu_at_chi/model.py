def run(inputs, settings):
    target = float(inputs["chi"])
    xs, ys = inputs["nu_curve"]["x"], inputs["nu_curve"]["y"]
    idx = min(range(len(xs)), key=lambda i: abs(xs[i] - target))
    if abs(xs[idx] - target) > 1e-4 * max(1.0, target):
        raise ValueError(f"chi={target} not present in curve x={xs}")
    return {"nu": float(ys[idx])}
