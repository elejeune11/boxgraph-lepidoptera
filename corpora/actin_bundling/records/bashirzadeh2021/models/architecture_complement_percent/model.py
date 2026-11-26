def run(inputs, settings):
    xs = inputs["ring_curve"]["x"]
    ys = inputs["ring_curve"]["y"]
    return {"non_ring_curve": {"x": [float(x) for x in xs], "y": [100.0 - float(y) for y in ys]}}
