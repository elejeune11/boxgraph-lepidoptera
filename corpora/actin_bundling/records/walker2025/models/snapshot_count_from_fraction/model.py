def run(inputs, settings):
    f_frac = float(inputs["f_last"]) / 100.0
    N_total = float(inputs["N_last"]) / f_frac
    dt_snap = float(inputs["t_sim"]) / N_total
    return {"N_total": N_total, "dt_snap": dt_snap}


def symbolic():
    import sympy as sp

    N_last, f_last, t_sim = sp.symbols("N_last f_last t_sim", positive=True)
    N_total = 100 * N_last / f_last
    dt_snap = t_sim / N_total
    return {"N_total": N_total, "dt_snap": dt_snap}
