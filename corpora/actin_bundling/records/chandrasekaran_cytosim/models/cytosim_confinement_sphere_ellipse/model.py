"""Faithful Python port of CytoSim spherical & ellipsoidal confinement:
- src/sim/spaces/space_sphere.cc:43-112 (SpaceSphere::volume, inside, setConfinement)
- src/sim/meca_inter.cc:4129, 4139-4167 (Meca::addSphereClamp)
- src/sim/spaces/space_ellipse.cc:89-191 (SpaceEllipse::volume, surface, inside, project3D)
- src/math/project_ellipse.cc:154-254 (projectEllipsoid Newton iteration)
- src/sim/space.cc:359-367 (Space::setConfinement -> Meca::addPlaneClamp)
- src/sim/meca_inter.cc:4764-4791 (Meca::addPlaneClamp)
"""
import math


def _project_ellipsoid(w, rad):
    # Line-by-line port of projectEllipsoid (src/math/project_ellipse.cc:154-254)
    aa = rad[0] * rad[0]
    bb = rad[1] * rad[1]
    cc = rad[2] * rad[2]

    h_min = 0.0
    if bb * cc * w[0] * w[0] + aa * cc * w[1] * w[1] + aa * bb * w[2] * w[2] < aa * bb * cc:
        h_min = math.sqrt(aa * w[0] * w[0] + bb * w[1] * w[1] + cc * w[2] * w[2]) - max(aa, bb, cc)

    h_min = max(h_min, abs(w[0]) * rad[0] - aa)
    h_min = max(h_min, abs(w[1]) * rad[1] - bb)
    h_min = max(h_min, abs(w[2]) * rad[2] - cc)

    h = h_min
    cnt = 0
    while True:
        aah = aa + h
        bbh = bb + h
        cch = cc + h

        wa_x = w[0] / aah
        wa_y = w[1] / bbh
        wa_z = w[2] / cch

        p_xx = wa_x * wa_x * aa
        p_yy = wa_y * wa_y * bb
        p_zz = wa_z * wa_z * cc

        h_old = h
        F = 1.0 - (p_xx + p_yy + p_zz)
        dF = 2.0 * (p_xx / aah + p_yy / bbh + p_zz / cch)
        if dF == 0.0:
            break
        h -= F / dF

        if h < h_min:
            h = 0.5 * (h_old + h_min)
            cnt += 1
            if cnt > 20:
                break
            continue

        cnt += 1
        if cnt > 20 or not (h > h_old):
            break

    return (
        w[0] * aa / (aa + h),
        w[1] * bb / (bb + h),
        w[2] * cc / (cc + h),
    )


def _ellipse_confinement_force_norm(w, rad, stiff: float) -> float:
    # Port of SpaceEllipse::inside (space_ellipse.cc:130) + Space::setConfinement (space.cc:359-367)
    inside = (w[0] / rad[0]) ** 2 + (w[1] / rad[1]) ** 2 + (w[2] / rad[2]) ** 2 <= 1.0
    if inside:
        return 0.0
    prj = _project_ellipsoid(w, rad)
    dx, dy, dz = prj[0] - w[0], prj[1] - w[1], prj[2] - w[2]
    n2 = dx * dx + dy * dy + dz * dz
    if n2 <= 1e-12:
        return 0.0
    # Meca::addPlaneClamp(mp, dir, -dot(prj, dir), stiff / n2):
    # At X = W, F = (stiff / n2) * dir * (dot(prj, dir) - dot(W, dir)) = stiff * (prj - W)
    return stiff * math.sqrt(n2)


def run(inputs, settings):
    R = float(inputs["R_sphere"])
    a = float(inputs["a_z"])
    b = float(inputs["b_xy"])
    exc_um = float(inputs["excursion"]) * 1e-3
    stiff = float(inputs["stiff"])

    # 1. SpaceSphere::setConfinement -> Meca::addSphereClamp (meca_inter.cc:4129)
    r_pos = R + exc_um
    if r_pos > R:
        # force = stiff * (C - P) * (1 - R / |PC|), norm = stiff * (r_pos - R)
        f_sphere = stiff * r_pos * (1.0 - R / r_pos)
    else:
        f_sphere = 0.0

    # 2. SpaceEllipse::setConfinement (inherited Space::setConfinement + projectEllipsoid)
    rad = (b, b, a)
    f_pole = _ellipse_confinement_force_norm((0.0, 0.0, a + exc_um), rad, stiff)
    f_equator = _ellipse_confinement_force_norm((b + exc_um, 0.0, 0.0), rad, stiff)
    # Off-axis exterior point (0.8*b, 0, 0.8*a): (0.8)^2 + (0.8)^2 = 1.28 > 1
    f_offaxis = _ellipse_confinement_force_norm((0.8 * b, 0.0, 0.8 * a), rad, stiff)

    vol_sphere = (4.0 * math.pi / 3.0) * (R ** 3)
    vol_ellipse = (4.0 * math.pi / 3.0) * b * b * a

    # Port of SpaceEllipse::surface() Knud Thomsen (space_ellipse.cc:112-118)
    POW = 1.6075
    ab = b * b
    ac = b * a
    bc = b * a
    s_thomsen = (4.0 * math.pi) * (((ab ** POW + ac ** POW + bc ** POW) / 3.0) ** (1.0 / POW))

    # Exact analytical prolate spheroid surface area (for comparison)
    if a > b:
        e = math.sqrt(1.0 - (b * b) / (a * a))
        s_exact = 2.0 * math.pi * b * b * (1.0 + (a / (b * e)) * math.asin(e))
    elif b > a:
        e = math.sqrt(1.0 - (a * a) / (b * b))
        s_exact = 2.0 * math.pi * b * b + math.pi * (a * a / e) * math.log((1.0 + e) / (1.0 - e))
    else:
        s_exact = 4.0 * math.pi * a * a

    return {
        "F_sphere": f_sphere,
        "F_ellipse_pole": f_pole,
        "F_ellipse_equator": f_equator,
        "F_ellipse_offaxis": f_offaxis,
        "volume_sphere": vol_sphere,
        "volume_ellipse": vol_ellipse,
        "surface_thomsen": s_thomsen,
        "surface_exact": s_exact,
        "aspect_ratio": a / b,
    }
