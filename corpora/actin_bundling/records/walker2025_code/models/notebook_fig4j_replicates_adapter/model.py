"""Live BoxGraph software adapter executing Figures.ipynb Cell 4 FO assembly unmodified."""

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
            next(reader)
            rows = [np.array([float(x) for x in r[1:]], dtype=float) for r in reader]
        return _MiniDF(rows)

    pd_mod.read_csv = read_csv
    return pd_mod


def run(inputs, settings):
    src_root = Path(os.environ["BOXGRAPH_SOURCE"])
    r_val = int(round(float(inputs["r_index"])))

    nb = json.loads((src_root / "Figures.ipynb").read_text(encoding="utf-8"))
    cell4_src = "".join(nb["cells"][4]["source"])
    ns = {"np": np, "pd": _make_pd_shim(src_root)}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(cell4_src, "Figures.ipynb#cell4", "exec"), ns)

    # Cell 4 builds a 3x3 grid over kbind in {0.1, 1.0, 10.0} (koniter=0..2 -> R row 2..4)
    # and kunbind in {0.01, 0.1, 1.0} (koffiter=0..2 -> R col 1..3):
    # Rval = (koniter + 2)*5 + (koffiter + 1), yielding local plot index idx = koniter*3 + koffiter.
    koniter = (r_val // 5) - 2
    koffiter = (r_val % 5) - 1
    if not (0 <= koniter < 3 and 0 <= koffiter < 3):
        raise ValueError(f"r_index={r_val} is outside Figure 4J's 3x3 subset")
    idx_first5 = koniter * 3 + koffiter
    idx_second5 = 9 + idx_first5

    # Cell 8 plots Meanlist[plotcounter*3 + val] for plotcounter, val in 0..2 -> indices 0..8 ONLY!
    meanlist = ns["Meanlist"]
    fo_first5 = float(meanlist[idx_first5][-1])
    fo_second5 = float(meanlist[idx_second5][-1])
    fo_all10 = 0.5 * (fo_first5 + fo_second5)

    # Number of replicates actually averaged in Meanlist[idx_first5] (5, since diriter=0 and 1 are separate entries)
    reps_per_entry = 5 if len(meanlist) == 18 else 10

    return {
        "replicates_plotted": float(reps_per_entry),
        "fo_t600_plotted": fo_first5,
        "fo_t600_second5": fo_second5,
        "fo_t600_all10": fo_all10,
    }
