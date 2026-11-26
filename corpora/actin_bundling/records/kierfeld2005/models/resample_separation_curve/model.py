import numpy as np


def run(inputs, settings):
    sx = np.array(inputs["source"]["x"], dtype=float)
    sy = np.array(inputs["source"]["y"], dtype=float)
    order = np.argsort(sx)
    sx, sy = sx[order], sy[order]
    lsx, lsy = np.log10(sx), np.log10(sy)
    tx = inputs["target"]["x"]
    x_min, x_max = float(sx[0]), float(sx[-1])
    ys_out = []
    for x in tx:
        if x_min * (1.0 - 1e-6) <= x <= x_max * (1.0 + 1e-6):
            val = float(10.0 ** np.interp(np.log10(x), lsx, lsy))
            ys_out.append(val)
        else:
            ys_out.append(None)
    return {"resampled": {"x": [float(x) for x in tx], "y": ys_out}}
