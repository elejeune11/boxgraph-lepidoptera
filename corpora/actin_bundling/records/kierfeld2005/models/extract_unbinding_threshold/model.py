def run(inputs, settings):
    xs = inputs["energy_curve"]["x"]
    ys = inputs["energy_curve"]["y"]
    pairs = sorted(zip(xs, ys), key=lambda p: p[0])
    # Find the discontinuous jump from unbound plateau (y > 0) to bound branch (y < 0)
    for (x1, y1), (x2, y2) in zip(pairs[:-1], pairs[1:]):
        if y1 > 0.0 and y2 < 0.0:
            return {"W_bar_c": 0.5 * (x1 + x2)}
    raise ValueError("No discontinuous drop across <H>/(NL) = 0 found in curve")
