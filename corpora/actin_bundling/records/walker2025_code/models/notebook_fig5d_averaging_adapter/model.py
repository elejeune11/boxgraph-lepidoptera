"""Live BoxGraph software adapter executing Figures.ipynb Cell 13 crossplothelper() unmodified."""

import contextlib
import csv
import io
import json
import os
import types
from pathlib import Path

import numpy as np


class _MiniIloc:
    def __init__(self, rows):
        self._rows = rows

    def __getitem__(self, idx):
        return _MiniSeries(self._rows[idx])


class _MiniSeries:
    def __init__(self, arr):
        self._arr = arr

    def to_numpy(self):
        return self._arr.copy()


class _MiniDF:
    def __init__(self, rows):
        self.iloc = _MiniIloc(rows)


def _make_pd_shim(src_root):
    pd_mod = types.ModuleType("pandas")

    def read_csv(filepath, index_col=0):
        full = Path(filepath)
        if not full.is_absolute():
            full = src_root / full
        with open(full, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            next(reader)  # header row
            rows = [np.array([float(x) for x in r[1:]], dtype=float) for r in reader]
        return _MiniDF(rows)

    pd_mod.read_csv = read_csv
    return pd_mod


def _extract_crossplothelper_source(nb_path):
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    cell_src = "".join(nb["cells"][13]["source"])
    # Extract only the crossplothelper definition (before the top-level loop separator)
    return cell_src.split("####")[0]


def run(inputs, settings):
    src_root = Path(os.environ["BOXGRAPH_SOURCE"])
    r_val = int(round(float(inputs["r_index"])))
    n_reps = int(round(float(inputs["n_reps"])))
    tv = float(inputs["n_monomers_total"])

    func_src = _extract_crossplothelper_source(src_root / "Figures.ipynb")
    observed_lens = []

    class _NpProxy:
        def __getattr__(self, name):
            if name == "mean":
                def _tracked_mean(a, *args, **kwargs):
                    observed_lens.append(int(np.asarray(a).size))
                    return np.mean(a, *args, **kwargs)
                return _tracked_mean
            return getattr(np, name)

    ns = {
        "np": _NpProxy(),
        "pd": _make_pd_shim(src_root),
        "Rval": r_val,
    }
    exec(compile(func_src, "Figures.ipynb#cell13", "exec"), ns)

    filename = f"outputfiles/CrosslinkOutput/Crosslink_inputfiles_R_{r_val}_r_"
    with contextlib.redirect_stdout(io.StringIO()):
        out = ns["crossplothelper"](filename, n_reps, tv)

    n_snaps_per_rep = float(observed_lens[0]) / float(n_reps)
    return {
        "n_snapshots_per_rep": n_snaps_per_rep,
        "frac_M0": float(out[4]),
        "frac_M1": float(out[5]),
        "frac_D0": float(out[6]),
        "frac_D1": float(out[7]),
        "frac_D2": float(out[8]),
    }
