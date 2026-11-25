import numpy as np


def run(inputs, settings):
    tau = np.asarray(inputs["DeltaL_vs_tau"]["x"], dtype=float)
    dL = np.asarray(inputs["DeltaL_vs_tau"]["y"], dtype=float)
    slope = float(np.polyfit(tau, dL, 1)[0])
    # 1 Pa*s * 1 um/s = 1 N/m^2 * s * 1e-6 m/s = 1e-6 N/m = 1 uN/m
    gamma = float(inputs["eta"]) * slope
    return {"slope": slope, "gamma": gamma}
