import numpy as np


def _peak_in_window(xs: np.ndarray, ys: np.ndarray, q_lo: float, q_hi: float) -> float:
    mask = (xs >= q_lo) & (xs <= q_hi)
    idx = np.where(mask)[0]
    if len(idx) == 0:
        raise ValueError(f"No curve samples in [{q_lo}, {q_hi}]")
    # Prefer true interior local maxima (avoids boundary shoulders of neighboring peaks)
    local_max = [int(i) for i in idx if 0 < i < len(xs) - 1 and ys[i] >= ys[i - 1] and ys[i] >= ys[i + 1]
                 and (ys[i] > ys[i - 1] or ys[i] > ys[i + 1])]
    cand = local_max if local_max else list(map(int, idx))
    best_i = max(cand, key=lambda i: float(ys[i]))
    # Local 3-point quadratic refinement when interior to the sampled curve and concave down
    if 0 < best_i < len(xs) - 1:
        x0, x1, x2 = float(xs[best_i - 1]), float(xs[best_i]), float(xs[best_i + 1])
        y0, y1, y2 = float(ys[best_i - 1]), float(ys[best_i]), float(ys[best_i + 1])
        if y1 >= y0 and y1 >= y2 and (y1 > y0 or y1 > y2):
            poly = np.polyfit([x0, x1, x2], [y0, y1, y2], 2)
            if poly[0] < 0.0:
                q_v = -poly[1] / (2.0 * poly[0])
                if x0 <= q_v <= x2:
                    return float(q_v)
    return float(xs[best_i])


def run(inputs, settings):
    curve = inputs["I_q"]
    xs = np.asarray(curve["x"], dtype=float)
    ys = np.asarray(curve["y"], dtype=float)
    order = np.argsort(xs)
    xs, ys = xs[order], ys[order]
    return {
        "q_layer1": _peak_in_window(xs, ys, 0.110, 0.124),
        "q_layer2": _peak_in_window(xs, ys, 0.124, 0.138),
    }
