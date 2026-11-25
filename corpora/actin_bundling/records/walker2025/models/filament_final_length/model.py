def run(inputs, settings):
    grown = float(inputs["L0"]) + float(inputs["k_grow"]) * float(inputs["t_sim"])
    return {"L_final": min(grown, float(inputs["L_max"]))}
