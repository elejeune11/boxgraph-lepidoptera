def run(inputs, settings):
    protein_ratio = inputs["I_spoke_protein"] / inputs["I_core_protein"]
    actin_ratio = inputs["I_spoke_actin"] / inputs["I_core_actin"]
    return {"enrichment": protein_ratio / actin_ratio}


def symbolic():
    import sympy as sp

    I_sp, I_cp, I_sa, I_ca = sp.symbols(
        "I_spoke_protein I_core_protein I_spoke_actin I_core_actin", positive=True
    )
    return {"enrichment": (I_sp / I_cp) / (I_sa / I_ca)}
