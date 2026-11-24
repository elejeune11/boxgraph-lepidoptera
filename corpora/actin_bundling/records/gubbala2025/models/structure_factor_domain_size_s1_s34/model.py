import numpy as np


def run(inputs, settings):
    k_max = float(inputs["k_max"])
    k_cut_ratio = float(inputs["k_cut_ratio"])
    L = float(inputs["L"])
    N = int(round(float(inputs["N"])))
    n_quad = int(settings.get("n_quad", 2000))
    rel_width = float(settings.get("rel_width", 0.25))

    k_cut = k_cut_ratio * k_max
    dk = 2.0 * np.pi / L
    # Evaluate Eq. (S1) & Eq. (S34) on a 2D-normalized test spectrum peaked at k_max
    # satisfying 2D Parseval normalization (1/(2*pi)) int k S(k) dk = 1 (equivalent to Eq. S33)
    k = np.linspace(dk, k_cut, n_quad)
    sigma = rel_width * k_max
    raw = np.exp(-0.5 * ((k - k_max) / sigma) ** 2)
    norm_2d = np.trapezoid(k * raw, k) / (2.0 * np.pi)
    S_k = raw / norm_2d if norm_2d > 0 else raw

    num = float(np.trapezoid(S_k, k))
    den = float(np.trapezoid(k * S_k, k))
    R = num / den
    k_m = den / num
    F_peak = float(k_m * k_m * np.max(S_k))
    return {"k_cut": k_cut, "k_m": k_m, "R": R, "F_peak": F_peak}
