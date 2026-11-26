import numpy as np


def _hex_coords(N):
    if N == 2:
        return np.array([[0.0, 0.0], [1.0, 0.0]])
    if N == 3:
        return np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3.0) / 2.0]])
    if N == 5:
        return np.array([
            [0.0, 0.0], [1.0, 0.0], [2.0, 0.0],
            [0.5, np.sqrt(3.0) / 2.0], [1.5, np.sqrt(3.0) / 2.0],
        ])
    pts = []
    lim = int(np.ceil(np.sqrt(N))) + 2
    for i in range(-lim, lim + 1):
        for j in range(-lim, lim + 1):
            x = i + 0.5 * j
            y = j * np.sqrt(3.0) / 2.0
            pts.append((x * x + y * y, x, y))
    pts.sort()
    return np.array([[p[1], p[2]] for p in pts[:N]])


def run(inputs, settings):
    N = int(round(inputs["N"]))
    l_r = inputs["l_r"]
    c = _hex_coords(N) * l_r
    drs = []
    for i in range(N):
        for j in range(i + 1, N):
            drs.append(float(np.linalg.norm(c[i] - c[j])) - l_r)
    return {"mean_separation": float(np.mean(drs))}
