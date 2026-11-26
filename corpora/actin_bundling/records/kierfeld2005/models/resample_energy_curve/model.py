import numpy as np


def run(inputs, settings):
    sx = np.array(inputs["source"]["x"], dtype=float)
    sy = np.array(inputs["source"]["y"], dtype=float)
    order = np.argsort(sx)
    sx, sy = sx[order], sy[order]
    tx = inputs["target"]["x"]
    x_min, x_max = float(sx[0]), float(sx[-1])
    ys_out = []
    for x in tx:
        if x_min - 1e-9 <= x <= x_max + 1e-9:
            ys_out.append(float(np.interp(x, sx, sy)))
        else:
            ys_out.append(None)
    return {"resampled": {"x": [float(x) for x in tx], "y": ys_out}}
