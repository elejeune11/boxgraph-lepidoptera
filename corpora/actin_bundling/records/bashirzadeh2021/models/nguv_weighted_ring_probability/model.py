def run(inputs, settings):
    r1 = inputs["ring_01"]["y"]
    r2 = inputs["ring_02"]["y"]
    r3 = inputs["ring_03"]["y"]
    xs = [float(x) for x in inputs["ring_01"]["x"]]
    n1 = [float(inputs["N01_1"]), float(inputs["N01_2"]), float(inputs["N01_3"])]
    n2 = [float(inputs["N02_1"]), float(inputs["N02_2"]), float(inputs["N02_3"])]
    n3 = [float(inputs["N03_1"]), float(inputs["N03_2"]), float(inputs["N03_3"])]
    out_y = []
    for i in range(3):
        num = n1[i] * float(r1[i]) + n2[i] * float(r2[i]) + n3[i] * float(r3[i])
        den = 100.0 * (n1[i] + n2[i] + n3[i])
        out_y.append(num / den)
    return {"pooled_ring_prob": {"x": xs, "y": out_y}}
