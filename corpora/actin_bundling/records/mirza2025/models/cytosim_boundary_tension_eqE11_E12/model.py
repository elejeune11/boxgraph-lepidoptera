def run(inputs, settings):
    L = inputs["L"]
    s_par = (inputs["F_top"] - inputs["F_bot"]) / (2.0 * L)
    s_perp = (inputs["F_right"] - inputs["F_left"]) / (2.0 * L)
    s_bar = 0.5 * (s_par + s_perp)
    xi = 2.0 * (s_par - s_perp) / (s_par + s_perp)
    return {
        "sigma_par": s_par,
        "sigma_perp": s_perp,
        "sigma_bar": s_bar,
        "xi": xi,
    }

def symbolic():
    import sympy as sp
    F_top, F_bot, F_right, F_left, L = sp.symbols("F_top F_bot F_right F_left L", real=True)
    s_par = (F_top - F_bot) / (2 * L)
    s_perp = (F_right - F_left) / (2 * L)
    s_bar = (s_par + s_perp) / 2
    xi = 2 * (s_par - s_perp) / (s_par + s_perp)
    return {"sigma_par": s_par, "sigma_perp": s_perp, "sigma_bar": s_bar, "xi": xi}
