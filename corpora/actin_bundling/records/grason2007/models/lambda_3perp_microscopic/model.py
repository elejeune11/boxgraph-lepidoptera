import math


def run(inputs, settings):
    # l_p is in um -> multiply by 1000 to convert to nm
    return {
        "lambda_3perp": math.sqrt(
            inputs["d"] * (1000.0 * inputs["l_p"]) / inputs["dE_kBT"]
        )
    }


def symbolic():
    import sympy as sp
    d, l_p, dE_kBT = sp.symbols("d l_p dE_kBT", positive=True)
    return {"lambda_3perp": sp.sqrt(d * 1000 * l_p / dE_kBT)}
