def run(inputs, settings):
    reg = inputs["I_region"]
    ext = inputs["I_exterior"]
    ext_by_x = {round(float(x), 8): float(y) for x, y in zip(ext["x"], ext["y"])}
    xs, ys = [], []
    for x, y in zip(reg["x"], reg["y"]):
        k = round(float(x), 8)
        if k in ext_by_x and ext_by_x[k] > 0:
            xs.append(float(x))
            ys.append(float(y) / ext_by_x[k])
    return {"RAM_aux": {"x": xs, "y": ys}}
