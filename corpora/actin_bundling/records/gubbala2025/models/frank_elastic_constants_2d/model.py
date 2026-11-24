def run(inputs, settings):
    E1 = float(inputs["E1"])
    E3 = float(inputs["E3"])
    S = float(inputs["S"])
    k11 = 2.0 * E1 * S**2 - E3 * S**3
    k33 = 2.0 * E1 * S**2 + E3 * S**3
    ratio = k33 / k11 if k11 != 0.0 else float("nan")
    return {"k11": k11, "k33": k33, "k33_over_k11": ratio}


def symbolic():
    import sympy as sp

    E1, E3, S = sp.symbols("E1 E3 S", real=True)
    k11 = 2 * E1 * S**2 - E3 * S**3
    k33 = 2 * E1 * S**2 + E3 * S**3
    return {"k11": k11, "k33": k33, "k33_over_k11": sp.simplify(k33 / k11)}
