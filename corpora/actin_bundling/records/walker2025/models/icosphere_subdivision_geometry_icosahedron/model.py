import math
import numpy as np


def _unit_icosahedron():
    t = (1.0 + math.sqrt(5.0)) / 2.0
    verts = np.array(
        [
            [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
            [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
            [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
        ],
        dtype=np.float64,
    )
    verts /= np.linalg.norm(verts, axis=1, keepdims=True)
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]
    return verts.tolist(), faces


def run(inputs, settings):
    r_nm = float(inputs["R_drop"]) * 1000.0
    n = int(round(float(inputs["n_subdiv"])))
    d_mon = float(inputs["d_mon"])

    verts, faces = _unit_icosahedron()
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
