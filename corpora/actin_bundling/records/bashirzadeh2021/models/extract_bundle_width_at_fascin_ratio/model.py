def run(inputs, settings):
    target = float(settings["target_ratio"])
    xs, ys = inputs["width_curve"]["x"], inputs["width_curve"]["y"]
    idx = min(range(len(xs)), key=lambda i: abs(xs[i] - target))
    if abs(xs[idx] - target) > 1e-4:
        raise ValueError(f"target_ratio={target} not present in curve x={xs}")
    return {"matched_width": {"x": [float(xs[idx])], "y": [float(ys[idx])]}}
