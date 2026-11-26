def run(inputs, settings):
    c_on = inputs["two_kon_over_D"]
    c_off = inputs["koff_over_D"]
    xs = [float(x) for x in c_on["x"]]
    ys = []
    for two_kon, koff in zip(c_on["y"], c_off["y"]):
        kon = float(two_kon) / 2.0
        ys.append((kon / float(koff)) ** 2)
    return {"N_peak": {"x": xs, "y": ys}}
