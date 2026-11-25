import numpy as np


def run(inputs, settings):
    prof = inputs["profile_C"]
    xs = np.array(prof["x"], dtype=float)
    ys = np.array(prof["y"], dtype=float)
    order = np.argsort(xs)
    xs, ys = xs[order], ys[order]
    c_edge = 0.5 * (float(ys[0]) + float(ys[-1]))
    c_half = 0.5 * (c_edge + float(np.max(ys)))
    above = np.where(ys >= c_half)[0]
    if len(above) < 2:
        return {"W": 0.0}
    iL, iR = int(above[0]), int(above[-1])
    if iL > 0 and ys[iL] != ys[iL - 1]:
        xL = xs[iL - 1] + (c_half - ys[iL - 1]) / (ys[iL] - ys[iL - 1]) * (xs[iL] - xs[iL - 1])
    else:
        xL = xs[iL]
    if iR < len(xs) - 1 and ys[iR + 1] != ys[iR]:
        xR = xs[iR] + (c_half - ys[iR]) / (ys[iR + 1] - ys[iR]) * (xs[iR + 1] - xs[iR])
    else:
        xR = xs[iR]
    return {"W": float(xR - xL)}
