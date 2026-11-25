"""Live BoxGraph software adapter calling Main_Analyze_Trajectory.find_chains_rings_and_doubles()."""

import importlib.util
import os
import sys
import types
from pathlib import Path

_MAT = None


def _load_main_analyze():
    global _MAT
    if _MAT is not None:
        return _MAT
    scripts_dir = Path(os.environ["BOXGRAPH_SOURCE"]) / "Analysis_Scripts"
    added_path = str(scripts_dir) not in sys.path
    if added_path:
        sys.path.insert(0, str(scripts_dir))
    injected_pd = "pandas" not in sys.modules
    if injected_pd:
        try:
            import pandas  # noqa: F401
            injected_pd = False
        except ImportError:
            sys.modules["pandas"] = types.ModuleType("pandas")
    try:
        spec = importlib.util.spec_from_file_location(
            "walker2025_main_analyze_topo", scripts_dir / "Main_Analyze_Trajectory.py"
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _MAT = mod
    finally:
        if injected_pd:
            sys.modules.pop("pandas", None)
        if added_path and str(scripts_dir) in sys.path:
            sys.path.remove(str(scripts_dir))
    return _MAT


def run(inputs, settings):
    mat = _load_main_analyze()
    chain_len = int(round(float(inputs["chain_len"])))
    ring_len = int(round(float(inputs["ring_len"])))
    n_doubles = int(round(float(inputs["n_doubles"])))

    mapping = {}
    next_id = 1

    # Linear chain of length chain_len
    if chain_len >= 2:
        ids = list(range(next_id, next_id + chain_len))
        next_id += chain_len
        for idx, nid in enumerate(ids):
            nbrs = []
            if idx > 0:
                nbrs.append(ids[idx - 1])
            if idx + 1 < chain_len:
                nbrs.append(ids[idx + 1])
            mapping[nid] = nbrs

    # Closed ring of length ring_len (>= 3)
    if ring_len >= 3:
        ids = list(range(next_id, next_id + ring_len))
        next_id += ring_len
        for idx, nid in enumerate(ids):
            mapping[nid] = [ids[(idx - 1) % ring_len], ids[(idx + 1) % ring_len]]

    # Double-bonded dimers (2 bonds between the same pair)
    for _ in range(n_doubles):
        a, b = next_id, next_id + 1
        next_id += 2
        mapping[a] = [b, b]
        mapping[b] = [a, a]

    chains, rings, doubles = mat.find_chains_rings_and_doubles(mapping)
    total_monomers = (
        sum(len(c) for c in chains)
        + sum(len(set(r)) for r in rings)
        + sum(len(d) for d in doubles)
    )
    return {
        "chains_found": float(len(chains)),
        "rings_found": float(len(rings)),
        "doubles_found": float(len(doubles)),
        "total_multimer_monomers": float(total_monomers),
    }
