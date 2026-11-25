from scipy import stats


def run(inputs, settings):
    t_val, p_val = stats.ttest_ind_from_stats(
        mean1=inputs["mean_1"],
        std1=inputs["sd_1"],
        nobs1=inputs["n_1"],
        mean2=inputs["mean_2"],
        std2=inputs["sd_2"],
        nobs2=inputs["n_2"],
        equal_var=True,
    )
    return {"t_stat": float(t_val), "p_value": float(p_val)}
