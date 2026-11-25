"""Live BoxGraph software adapter auditing Main_Analyze_Trajectory.getvaspprops() in memory."""

import contextlib
import importlib.util
import io
import os
import sys
import types
from pathlib import Path

import numpy as np


def _make_pandas_capture(captured):
    pd_mod = types.ModuleType("pandas")

    class _FakeDF:
        def __init__(self, data, index=None):
            self.data = np.array(data, dtype=float)
            self.index = list(index) if index is not None else []

        def to_csv(self, path):
            captured["data"] = self.data
            captured["index"] = self.index

    pd_mod.DataFrame = _FakeDF
    return pd_mod


def run(inputs, settings):
    n_cross = int(round(float(inputs["n_crosslinkers"])))
    n_snaps = int(settings.get("n_snaps", 20))
    scripts_dir = Path(os.environ["BOXGRAPH_SOURCE"]) / "Analysis_Scripts"
    captured = {}

    added_path = str(scripts_dir) not in sys.path
    if added_path:
        sys.path.insert(0, str(scripts_dir))
    orig_pd = sys.modules.get("pandas")
    sys.modules["pandas"] = _make_pandas_capture(captured)
    try:
        spec = importlib.util.spec_from_file_location(
            "walker2025_main_analyze_vasp", scripts_dir / "Main_Analyze_Trajectory.py"
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        # Build an in-memory replicate where ALL n_cross crosslinkers bridge 2 distinct filaments (1 & 2)
        rep = mod.replicate(0)
        for s in range(n_snaps):
            rep.addsnapshot()
            snap = rep.snap[s]
            for cid in range(1, n_cross + 1):
                snap.crossboundblobid.extend([cid, cid])
                snap.crossboundfilid.extend([1, 2])

        with contextlib.redirect_stdout(io.StringIO()):
            mod.getvaspprops([rep], 1, n_cross, 0, "in_memory_only.csv")
    finally:
        if orig_pd is None:
            sys.modules.pop("pandas", None)
        else:
            sys.modules["pandas"] = orig_pd
        if added_path and str(scripts_dir) in sys.path:
            sys.path.remove(str(scripts_dir))

    idx = captured["index"]
    mat = captured["data"]
    nfpc_2 = float(np.mean(mat[idx.index("Mean_Nfpc_2_fil")]))
    nfpc_4 = float(np.mean(mat[idx.index("Mean_Nfpc_4_fil")]))
    return {
        "mean_nfpc_2_fil": nfpc_2,
        "mean_nfpc_4_fil": nfpc_4,
        "frac_misclassified_to_4fil": nfpc_4 / float(n_cross),
    }
