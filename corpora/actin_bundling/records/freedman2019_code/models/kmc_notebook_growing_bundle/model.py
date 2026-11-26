"""Live BoxGraph software adapter calling extension_population_ratio() directly from xlink_sorting_paper/kmc.ipynb."""

import contextlib
import io
import json
import math
import os
from pathlib import Path
import numpy as np


def _load_kmc_namespace():
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "xlink_sorting_paper" / "kmc.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    ns = {"np": np, "exp": np.exp, "log": np.log, "sqrt": np.sqrt}
    # Cell 3: f() and ratio_at_dist()
    src3 = "".join(nb["cells"][3]["source"])
    exec(compile(src3, "kmc.ipynb#cell3", "exec"), ns)
    # Cell 5: rates() and extension_population_ratio(), omitting the trailing interactive demo call
    lines5 = [ln for ln in nb["cells"][5]["source"] if not ln.startswith("extension_population_ratio(")]
    src5 = "".join(lines5)
    exec(compile(src5, "kmc.ipynb#cell5", "exec"), ns)
    return ns


def run(inputs, settings):
    ns = _load_kmc_namespace()
    site_width = float(inputs["site_width"])
    dxl = float(inputs["dxl"])
    max_len = float(inputs["max_filament_length"])
    ext_rate = float(inputs["extension_rate"])
    kon2 = float(inputs["kon2"])
    kratio = float(inputs["kratio"])
    kon1 = kon2 * kratio
    Lp = float(inputs["Lp"])

    ns["site_width"] = site_width
    ns["dxl"] = dxl

    seed = int(settings.get("seed", 1))
    num_trials = int(settings.get("num_trials", 1000))
    max_extensions = int(settings.get("max_extensions", 200))
    warmup_trial = bool(settings.get("warmup_trial", False))

    np.random.seed(seed)
    with contextlib.redirect_stdout(io.StringIO()):
        if warmup_trial:
            # Matches kmc.ipynb Cell 5's trailing num_trials=1 demo call prior to Cell 6
            ns["extension_population_ratio"](
                num_trials=1,
                max_filament_length=max_len,
                extension_rate=ext_rate,
                kon1=kon1,
                kon2=kon2,
                showplot=False,
                Lp=Lp,
                max_extensions=max_extensions,
            )
        swap_ratio, _, domain_lengths, domain_types, _, _, _, _ = ns["extension_population_ratio"](
            num_trials=num_trials,
            max_filament_length=max_len,
            extension_rate=ext_rate,
            kon1=kon1,
            kon2=kon2,
            showplot=False,
            Lp=Lp,
            max_extensions=max_extensions,
        )

    swap_ratio = float(swap_ratio)
    eff_barrier = -math.log(swap_ratio) if swap_ratio > 0.0 else float("inf")
    dl_a = float(domain_lengths[domain_types == 1].mean()) / 1000.0 if np.any(domain_types == 1) else 0.0
    dl_b = float(domain_lengths[domain_types == 2].mean()) / 1000.0 if np.any(domain_types == 2) else 0.0
    max_n_sites = int(max_len / site_width)

    return {
        "switching_ratio": swap_ratio,
        "effective_barrier": eff_barrier,
        "mean_domain_length_A": dl_a,
        "mean_domain_length_B": dl_b,
        "max_n_sites": float(max_n_sites),
    }
