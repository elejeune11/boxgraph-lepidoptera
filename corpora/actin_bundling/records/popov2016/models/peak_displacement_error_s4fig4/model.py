def run(inputs, settings):
    ys = inputs["err_curve"]["y"]
    return {"max_e_d": float(max(ys))}
