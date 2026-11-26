def run(inputs, settings):
    curve = inputs["D_b_mean"]
    D_f = float(inputs["D_f"])
    xs = [float(x) for x in curve["x"]]
    ys = [(float(d) / D_f) ** 2 for d in curve["y"]]
    return {"N_b_mean": {"x": xs, "y": ys}}
