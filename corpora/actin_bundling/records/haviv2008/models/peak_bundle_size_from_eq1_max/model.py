import math
from scipy.optimize import minimize_scalar


def run(inputs, settings):
    kon_D = float(inputs["kon_over_D"])
    koff_D = float(inputs["koff_over_D"])
    xatol = float(settings.get("xatol", 1.0e-7))

    # Bracket cleanly above the location where the derivative turns negative
    hi = max(10.0, 10.0 * (kon_D / koff_D) ** 2 + 10.0)

    def neg_log_P(N_b):
        return -(2.0 * kon_D * math.sqrt(N_b) - koff_D * N_b)

    res = minimize_scalar(
        neg_log_P,
        bounds=(1.0e-12, hi),
        method="bounded",
        options={"xatol": xatol},
    )
    return {"N_peak": float(res.x)}
