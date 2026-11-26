def run(inputs, settings):
    n = int(settings.get("n_profilin", 2))
    return {"MW_complex": inputs["MW_vasp"] + n * inputs["MW_profilin"]}
