"""Live BoxGraph software adapter calling f() and ratio_at_dist() from xlink_sorting_paper/kmc.ipynb."""

import json
import os
from pathlib import Path
import numpy as np


def _load_kmc_cells():
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "xlink_sorting_paper" / "kmc.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    ns = {"np": np, "exp": np.exp, "log": np.log, "sqrt": np.sqrt}
    # Load Cell 3 (f and ratio_at_dist definitions) directly from the released notebook
    src3 = "".join(nb["cells"][3]["source"])
    exec(compile(src3, "kmc.ipynb#cell3", "exec"), ns)
    return ns


def run(inputs, settings):
    ns = _load_kmc_cells()
    site_width = float(inputs["site_width"])
    dxl = float(inputs["l_long"]) - float(inputs["l_short"])
    Lp = float(inputs["Lp"])
    x = float(inputs["x_sites"])
    max_len = float(inputs["max_filament_length"])

    # Inject module-level parameters expected by kmc.ipynb's f() and ratio_at_dist()
    ns["site_width"] = site_width
    ns["dxl"] = dxl

    u_kbt = float(ns["f"](x, Lp))
    dist_ratio = float(ns["ratio_at_dist"](x, onleft=1, kon1=1.0, kon2=1.0, Lp=Lp))
    p_switch = dist_ratio / (1.0 + dist_ratio)
    max_n_sites = int(max_len / site_width)
    gap_distance_um = (x * site_width) / 1000.0

    return {
        "dxl": dxl,
        "max_n_sites": float(max_n_sites),
        "gap_distance": gap_distance_um,
        "U_kBT": u_kbt,
        "p_switch": p_switch,
    }
