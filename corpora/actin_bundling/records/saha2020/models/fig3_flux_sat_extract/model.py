import numpy as np


def run(inputs, settings):
    ts = np.array([float(v) for v in inputs["J_curve"]["x"]])
    Js = np.array([float(v) for v in inputs["J_curve"]["y"]])
    t_sat = float(settings.get("t_sat", 4.0))
    m = ts >= t_sat
    return {"J_sat": float(np.mean(Js[m]))}
