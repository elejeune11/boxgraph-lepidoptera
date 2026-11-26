"""Sum of mutually exclusive membrane-localized actin morphology percentages."""


def run(inputs, settings):
    return {"f_total": float(inputs["f_domains"]) + float(inputs["f_cortex"])}


def symbolic():
    import sympy as sp

    f_domains, f_cortex = sp.symbols("f_domains f_cortex", nonnegative=True)
    return {"f_total": f_domains + f_cortex}
