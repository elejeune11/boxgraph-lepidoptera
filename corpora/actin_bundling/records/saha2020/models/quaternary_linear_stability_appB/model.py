import math
import numpy as np


def run(inputs, settings):
    c = np.array([inputs["c1"], inputs["c2"], inputs["c3"], inputs["c4"]], dtype=float)
    pb = np.array([inputs["phi1_bar"], inputs["phi2_bar"], inputs["phi3_bar"], inputs["phi4_bar"]], dtype=float)
    kappa = float(inputs["kappa"])

    C = np.array([
        [0.0, inputs["C12"], inputs["C13"], inputs["C14"]],
        [inputs["C21"], 0.0, inputs["C23"], inputs["C24"]],
        [inputs["C31"], inputs["C32"], 0.0, inputs["C34"]],
        [inputs["C41"], inputs["C42"], inputs["C43"], 0.0],
    ], dtype=float)

    M = -C.copy()
    for i in range(4):
        M[i, i] = 4.0 * (c[i] ** 2 - 3.0 * pb[i] ** 2)

    eigs = np.linalg.eigvals(M)
    # Collapse conjugate pairs by keeping eigenvalues with Im >= -1e-12
    unique = []
    for e in eigs:
        if e.imag >= -1e-10:
            unique.append(e)
    unique.sort(key=lambda z: z.real, reverse=True)

    lead = unique[0]
    sec = unique[1] if len(unique) > 1 else unique[0]

    lam1 = 2.0 * math.pi / math.sqrt(lead.real / (2.0 * kappa)) if lead.real > 0 else float("nan")
    lam2 = 2.0 * math.pi / math.sqrt(sec.real / (2.0 * kappa)) if sec.real > 0 else float("nan")
    return {
        "re_leading_over_q2": float(lead.real),
        "im_leading_over_q2": float(abs(lead.imag)),
        "wavelength_leading": float(lam1),
        "wavelength_secondary": float(lam2),
    }
