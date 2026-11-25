def run(inputs, settings):
    ys = inputs["theta_vs_Rmax"]["y"]
    return {"theta_mean": sum(ys) / len(ys)}
