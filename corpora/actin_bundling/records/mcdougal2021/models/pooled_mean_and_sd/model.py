"""Pooled mean and pooled standard deviation across 3 scale profiles (p. 7, SI p. 11).

Computes both the exact combined sample SD of the pooled dataset (pooled_sd_total)
and the ANOVA within-group pooled SD (pooled_sd_within) to make explicit the
statistical convention not written out in the paper text.
"""

import math


def run(inputs, settings):
    ns = [float(inputs["n1"]), float(inputs["n2"]), float(inputs["n3"])]
    mus = [float(inputs["mu1"]), float(inputs["mu2"]), float(inputs["mu3"])]
    sds = [float(inputs["sd1"]), float(inputs["sd2"]), float(inputs["sd3"])]

    N = sum(ns)
    K = len(ns)
    mu_pool = sum(n * mu for n, mu in zip(ns, mus)) / N

    ss_within = sum((n - 1.0) * (s * s) for n, s in zip(ns, sds))
    ss_between = sum(n * (mu - mu_pool) ** 2 for n, mu in zip(ns, mus))

    s_total = math.sqrt((ss_within + ss_between) / (N - 1.0))
    s_within = math.sqrt(ss_within / (N - K))
    return {
        "pooled_mean": mu_pool,
        "pooled_sd_total": s_total,
        "pooled_sd_within": s_within,
    }


def symbolic():
    import sympy as sp

    n1, n2, n3 = sp.symbols("n1 n2 n3", positive=True)
    mu1, mu2, mu3 = sp.symbols("mu1 mu2 mu3", real=True)
    sd1, sd2, sd3 = sp.symbols("sd1 sd2 sd3", nonnegative=True)

    N = n1 + n2 + n3
    mu_pool = (n1 * mu1 + n2 * mu2 + n3 * mu3) / N
    ss_w = (n1 - 1) * sd1**2 + (n2 - 1) * sd2**2 + (n3 - 1) * sd3**2
    ss_b = n1 * (mu1 - mu_pool) ** 2 + n2 * (mu2 - mu_pool) ** 2 + n3 * (mu3 - mu_pool) ** 2
    return {
        "pooled_mean": mu_pool,
        "pooled_sd_total": sp.sqrt((ss_w + ss_b) / (N - 1)),
        "pooled_sd_within": sp.sqrt(ss_w / (N - 3)),
    }
