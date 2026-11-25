"""Upper-tail chi-squared p-value for likelihood ratio tests (p. 3, p. 13)."""

from scipy import stats


def run(inputs, settings):
    chi2 = float(inputs["chi2"])
    df = float(inputs["df"])
    return {"p_value": float(stats.chi2.sf(chi2, df))}
