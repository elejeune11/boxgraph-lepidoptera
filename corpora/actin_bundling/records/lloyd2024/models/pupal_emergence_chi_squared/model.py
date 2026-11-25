"""Pupal emergence accounting and 2x2 contingency table chi-squared test (p. 13, SI p. 6)."""

from scipy import stats


def run(inputs, settings):
    n_t = float(inputs["n_treated_injected"])
    n_c = float(inputs["n_control_injected"])
    e_t = float(inputs["n_treated_emerged"])
    e_c = float(inputs["n_control_emerged"])

    n_tot = n_t + n_c
    e_tot = e_t + e_c
    rate_pct = 100.0 * e_tot / n_tot

    table = [
        [e_t, n_t - e_t],
        [e_c, n_c - e_c],
    ]
    correction = bool(settings.get("yates_correction", False))
    chi2, p_val, _, _ = stats.chi2_contingency(table, correction=correction)

    return {
        "n_total_injected": n_tot,
        "emergence_rate": rate_pct,
        "chi2": float(chi2),
        "p_value": float(p_val),
    }
