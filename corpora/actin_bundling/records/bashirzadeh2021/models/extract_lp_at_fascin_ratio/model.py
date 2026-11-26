def run(inputs, settings):
    target = float(settings["target_ratio"])
    xs, ys = inputs["lp_curve"]["x"], inputs["lp_curve"]["y"]
    idx = min(range(len(xs)), key=lambda i: abs(xs[i] - target))
    if abs(xs[idx] - target) > 1e-4:
        raise ValueError(f"target_ratio={target} not present in curve x={xs}")
    return {"L_p": float(ys[idx])}
