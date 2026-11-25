def run(inputs, settings):
    ys = inputs["theta_vs_n"]["y"]
    return {"theta_mean": sum(ys) / len(ys)}
