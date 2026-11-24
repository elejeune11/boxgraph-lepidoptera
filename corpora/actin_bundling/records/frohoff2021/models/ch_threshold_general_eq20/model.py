def run(inputs, settings):
    f2_pp = inputs["f2_pp"]
    if abs(f2_pp) < 1e-15:
        return {"f1_pp_CH": None}
    Delta = inputs["alpha"] ** 2 - inputs["rho"] ** 2
    return {"f1_pp_CH": -Delta / f2_pp}


def symbolic():
    import sympy as sp

    f2_pp, rho, alpha = sp.symbols("f2_pp rho alpha", real=True)
    return {"f1_pp_CH": -(alpha**2 - rho**2) / f2_pp}
