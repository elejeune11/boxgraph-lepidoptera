def run(inputs, settings):
    xs, ys = inputs["length_curve"]["x"], inputs["length_curve"]["y"]
    idx = max(range(len(xs)), key=lambda i: xs[i])
    return {"L_final": float(ys[idx])}
