def run(inputs, settings):
    return {"var": inputs["mean_sq_disp"] - inputs["mean_abs_disp"] ** 2}


def symbolic():
    import sympy as sp

    mean_sq_disp, mean_abs_disp = sp.symbols("mean_sq_disp mean_abs_disp", nonnegative=True)
    return {"var": mean_sq_disp - mean_abs_disp**2}
