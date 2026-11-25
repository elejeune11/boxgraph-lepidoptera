def run(inputs, settings):
    return {"mc_total": inputs["mc_base"] + inputs["mc_added"]}


def symbolic():
    import sympy as sp
    mc_base, mc_added = sp.symbols("mc_base mc_added", positive=True)
    return {"mc_total": mc_base + mc_added}
