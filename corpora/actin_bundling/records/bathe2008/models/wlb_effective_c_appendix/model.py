"""Appendix Eqs. (A8)-(A10) effective mode factor c_eff = [sum c(q_j)/q_j^4] / [sum 1/q_j^4]."""

import math


def run(inputs, settings):
    bc = settings.get("bc", "pinned")
    n_modes = int(settings.get("n_modes", 2000))
    if bc == "pinned":
        qLs = [(2 * k + 1) * math.pi for k in range(n_modes)]
    elif bc == "clamped":
        qLs = [2.0 * (2 * k + 1) * math.pi for k in range(n_modes)]
    elif bc == "ring":
        qLs = [2.0 * (k + 1) * math.pi for k in range(n_modes)]
    else:
        raise ValueError(f"Unknown bc: {bc}")
    num = sum(((qL ** 2) / 12.0) / (qL ** 4) for qL in qLs)
    den = sum(1.0 / (qL ** 4) for qL in qLs)
    return {"c_eff": num / den}
