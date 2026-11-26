def run(inputs, settings):
    area_bin = (inputs["L_domain"] / inputs["N_bins"]) ** 2
    return {"pdf_val": inputs["n_bin"] / (area_bin * inputs["N_total"])}


def symbolic():
    import sympy as sp

    n_bin, N_total, L_domain, N_bins = sp.symbols("n_bin N_total L_domain N_bins", positive=True)
    return {"pdf_val": n_bin / ((L_domain / N_bins) ** 2 * N_total)}
