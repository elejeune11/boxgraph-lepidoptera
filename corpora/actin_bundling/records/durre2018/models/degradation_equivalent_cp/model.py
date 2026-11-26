"""Equivalent CP concentration matching the first-order barbed-end degradation rate k_Deg: k_Cap * CP_eq = k_Deg."""


def run(inputs, settings):
    k_deg = float(inputs["k_Deg"])  # 1/s
    k_cap = float(inputs["k_Cap"])  # 1/(uM*s)
    # uM -> nM factor is 1000
    return {"CP_eq": (k_deg / k_cap) * 1000.0}


def symbolic():
    import sympy as sp

    k_Deg, k_Cap = sp.symbols("k_Deg k_Cap", positive=True)
    return {"CP_eq": 1000 * k_Deg / k_Cap}
