def run(inputs, settings):
    k5 = float(inputs["k_5"])
    m5 = float(inputs["m_5"])
    n5 = float(inputs["n_5"])
    t5 = float(inputs["t_5"])
    k4 = float(inputs["k_4"])
    m4 = float(inputs["m_4"])
    n4 = float(inputs["n_4"])
    t4 = float(inputs["t_4"])
    N_repeat = k5 * m5 + k4 * m4
    T_repeat = k5 * m5 * (t5 / n5) + k4 * m4 * (t4 / n4)
    f_bound = (k5 + k4) / N_repeat
    return {
        "N_repeat": N_repeat,
        "T_repeat": T_repeat,
        "f_bound": f_bound,
    }


def symbolic():
    import sympy as sp

    k_5, m_5, n_5, t_5, k_4, m_4, n_4, t_4 = sp.symbols(
        "k_5 m_5 n_5 t_5 k_4 m_4 n_4 t_4", positive=True
    )
    N_repeat = k_5 * m_5 + k_4 * m_4
    T_repeat = k_5 * m_5 * t_5 / n_5 + k_4 * m_4 * t_4 / n_4
    return {
        "N_repeat": N_repeat,
        "T_repeat": T_repeat,
        "f_bound": (k_5 + k_4) / N_repeat,
    }
