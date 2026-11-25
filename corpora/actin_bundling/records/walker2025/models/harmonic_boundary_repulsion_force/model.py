def run(inputs, settings):
    overshoot = max(0.0, float(inputs["r"]) - float(inputs["R_drop"]))
    return {"F_rep": float(inputs["k_bound"]) * overshoot}
