def run(inputs, settings):
    return {
        "pct_total": (
            float(inputs["pct_DOPC"])
            + float(inputs["pct_NiNTA"])
            + float(inputs["pct_biotin"])
            + float(inputs["pct_TR"])
        )
    }


def symbolic():
    import sympy as sp

    pct_DOPC, pct_NiNTA, pct_biotin, pct_TR = sp.symbols(
        "pct_DOPC pct_NiNTA pct_biotin pct_TR", positive=True
    )
    return {"pct_total": pct_DOPC + pct_NiNTA + pct_biotin + pct_TR}
