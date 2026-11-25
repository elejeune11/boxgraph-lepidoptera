import math
import numpy as np


def _unit_octahedron():
    verts = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]
    faces = [
        (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
        (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5),
    ]
    return verts, faces


def run(inputs, settings):
    r_nm = float(inputs["R_drop"]) * 1000.0
    n = int(round(float(inputs["n_subdiv"])))
    d_mon = float(inputs["d_mon"])

    verts, faces = _unit_octahedron()
    for _ in range(n):
        cache = {}
        new_faces = []

        def get_mid(i, j):
            key = (i, j) if i < j else (j, i)
            if key in cache:
                return cache[key]
            v1, v2 = verts[i], verts[j]
            m = [0.5 * (v1[k] + v2[k]) for k in range(3)]
            norm = math.sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2])
            verts.append([m[0] / norm, m[1] / norm, m[2] / norm])
            idx = len(verts) - 1
            cache[key] = idx
            return idx

        for a, b, c in faces:
            ab = get_mid(a, b)
            bc = get_mid(b, c)
            ca = get_mid(c, a)
            new_faces.extend([(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)])
        faces = new_faces

    v_arr = np.array(verts, dtype=np.float64) * r_nm
    edges = set()
    for a, b, c in faces:
        edges.add((a, b) if a < b else (b, a))
        edges.add((b, c) if b < c else (c, b))
        edges.add((c, a) if c < a else (a, c))
    edge_idx = np.array(list(edges), dtype=np.int64)
    lengths = np.linalg.norm(v_arr[edge_idx[:, 0]] - v_arr[edge_idx[:, 1]], axis=1)

    n_faces = len(faces)
    face_area = (4.0 * math.pi * r_nm * r_nm) / n_faces
    eq_side = math.sqrt(4.0 * face_area / math.sqrt(3.0))

    return {
        "n_faces": float(n_faces),
        "median_edge": float(np.median(lengths)),
        "eq_side": eq_side,
        "n_monomers": eq_side / d_mon,
    }
