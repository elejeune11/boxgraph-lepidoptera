"""Welch two-sample t-test from group means, standard errors, and sample sizes (p. 10, SI p. 6)."""

import math
from scipy import stats


def run(inputs, settings):
    m1 = float(inputs["mean_control"])
    se1 = float(inputs["se_control"])
    n1 = float(inputs["n_control"])
    m2 = float(inputs["mean_treated"])
    se2 = float(inputs["se_treated"])
    n2 = float(inputs["n_treated"])

    v1 = se1 * se1
    v2 = se2 * se2
    denom = math.sqrt(v1 + v2)
    t_stat = (m1 - m2) / denom
    df = (v1 + v2) ** 2 / ((v1 ** 2) / (n1 - 1.0) + (v2 ** 2) / (n2 - 1.0))
    p_val = float(2.0 * stats.t.sf(abs(t_stat), df))
    return {"t_stat": t_stat, "df": df, "p_value": p_val}
