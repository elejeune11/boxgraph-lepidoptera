import math


def run(inputs, settings):
    eps = inputs["epsilon_kBT"]
    ratio = math.exp(-eps)
    return {
        "mu_on_diff": inputs["mu_on_same"] * ratio,
        "switching_ratio": ratio,
        "fold_reduction": math.exp(eps),
    }


def symbolic():
    import sympy as sp

    mu_on_same, epsilon_kBT = sp.symbols("mu_on_same epsilon_kBT", positive=True)
    ratio = sp.exp(-epsilon_kBT)
    return {
        "mu_on_diff": mu_on_same * ratio,
        "switching_ratio": ratio,
        "fold_reduction": sp.exp(epsilon_kBT),
    }
