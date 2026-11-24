def run(inputs, settings):
    return {"dm_dt": -inputs["k_off"] * inputs["m"] + inputs["k_on"] * inputs["n"]}


def symbolic():
    import sympy as sp

    k_off, k_on, m, n = sp.symbols("k_off k_on m n", real=True)
    return {"dm_dt": -k_off * m + k_on * n}
