import numpy as np


def run(inputs, settings):
    r_c = float(inputs["r_c"])
    sigma = float(inputs["sigma"])
    n_bins = int(round(settings.get("n_bins", 20)))

    edges = np.linspace(-r_c, r_c, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    X, Y = np.meshgrid(centers, centers)
    R2 = X ** 2 + Y ** 2
    area_bin = ((2.0 * r_c) / n_bins) ** 2

    raw = np.where(R2 <= r_c ** 2, np.exp(-R2 / (2.0 * sigma ** 2)), 0.0)
    s = raw / (np.sum(raw) * area_bin)

    s2_sum = float(np.sum(s ** 2))
    s4_sum = float(np.sum(s ** 4))
    U_4 = 1.0 - s4_sum / (3.0 * (s2_sum ** 2))
    return {"U_4": float(U_4), "peak_pdf": float(np.max(s))}
