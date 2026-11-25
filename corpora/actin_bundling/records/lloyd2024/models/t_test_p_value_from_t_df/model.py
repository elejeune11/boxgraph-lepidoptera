"""Two-sided Student t-distribution p-value from printed t statistic and df."""

from scipy import stats


def run(inputs, settings):
    t_stat = float(inputs["t_stat"])
    df = float(inputs["df"])
    return {"p_value": float(2.0 * stats.t.sf(abs(t_stat), df))}
