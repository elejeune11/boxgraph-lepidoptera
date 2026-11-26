def run(inputs, settings):
    xs, ys = inputs["conc_curve"]["x"], inputs["conc_curve"]["y"]
    idx = max(range(len(xs)), key=lambda i: xs[i])
    return {"A_final": float(ys[idx])}
