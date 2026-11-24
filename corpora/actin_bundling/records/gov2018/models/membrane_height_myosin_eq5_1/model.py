def run(inputs, settings):
    dh_dt = (
        inputs["F_passive"]
        + inputs["A"] * (inputs["n"] - inputs["n0"])
        - inputs["A_star"] * (inputs["m"] - inputs["m0"])
    ) / inputs["xi"]
    return {"dh_dt": dh_dt}


def symbolic():
    import sympy as sp

    xi, F_passive, A, n, n0, A_star, m, m0 = sp.symbols(
        "xi F_passive A n n0 A_star m m0", real=True
    )
    return {"dh_dt": (F_passive + A * (n - n0) - A_star * (m - m0)) / xi}
