import math
from scipy import stats


def run(inputs, settings):
    m1 = float(inputs["mean1"])
    se1 = float(inputs["sem1"])
    n1 = float(inputs["n1"])
    m2 = float(inputs["mean2"])
    se2 = float(inputs["sem2"])
    n2 = float(inputs["n2"])

    v1 = se1 ** 2
    v2 = se2 ** 2
    t_stat = (m2 - m1) / math.sqrt(v1 + v2)
    df = ((v1 + v2) ** 2) / ((v1 ** 2) / (n1 - 1.0) + (v2 ** 2) / (n2 - 1.0))
    p_val = float(2.0 * stats.t.sf(abs(t_stat), df))
    return {"p_value": p_val}
