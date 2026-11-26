def run(inputs, settings):
    return {"ratio": inputs["tau_m_long"] / inputs["tau_m_short"]}


def symbolic():
    import sympy as sp

    tau_m_long, tau_m_short = sp.symbols("tau_m_long tau_m_short", positive=True)
    return {"ratio": tau_m_long / tau_m_short}
