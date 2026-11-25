"""Live BoxGraph software adapter calling icosphere.icosphere() directly from $BOXGRAPH_SOURCE."""

import importlib.util
import math
import os
from pathlib import Path

import numpy as np

_ICO = None


def _load_icosphere():
    global _ICO
    if _ICO is not None:
        return _ICO
    src = Path(os.environ["BOXGRAPH_SOURCE"]) / "Analysis_Scripts" / "icosphere.py"
    spec = importlib.util.spec_from_file_location("walker2025_icosphere", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _ICO = mod
    return _ICO


def run(inputs, settings):
    ic = _load_icosphere()
    R_um = float(inputs["R_drop"])
    subdiv = int(round(float(inputs["n_subdiv"])))
    d_mon_nm = float(inputs["d_mon"])

    verts, faces = ic.icosphere(subdiv, R_um)
    n_faces = int(faces.shape[0])
    n_vertices = int(verts.shape[0])

    # Chordal edge lengths across all mesh triangles (um -> nm)
    v0 = verts[faces[:, 0]]
    v1 = verts[faces[:, 1]]
    v2 = verts[faces[:, 2]]
    e01 = np.linalg.norm(v1 - v0, axis=1)
    e12 = np.linalg.norm(v2 - v1, axis=1)
    e20 = np.linalg.norm(v0 - v2, axis=1)
    median_edge_nm = float(np.median(np.concatenate([e01, e12, e20])) * 1000.0)

    # Planar triangle mean area -> area-equivalent equilateral side (nm)
    tri_areas_um2 = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    mean_area_nm2 = float(np.mean(tri_areas_um2) * 1e6)
    eq_side_nm = math.sqrt(4.0 * mean_area_nm2 / math.sqrt(3.0))
    n_monomers = eq_side_nm / d_mon_nm

    return {
        "n_faces": float(n_faces),
        "n_vertices": float(n_vertices),
        "median_edge": float(median_edge_nm),
        "eq_side": float(eq_side_nm),
        "n_monomers": float(n_monomers),
    }
