"""Faithful Python port of filament::grow() and filament::update_length() in src/filament.cpp:673-736.

Simulates the exact floating-point split rule:
    if (lb + dL < l0_max): lb = lb + dL
    else: split barbed-end link into two links of rest length spring_l0
and computes the exact realized contour growth rate and discarded overshoot per split.
"""


def run(inputs, settings):
    kgrow = float(inputs["kgrow"])
    dL = float(inputs["lgrow"])
    spring_l0 = float(inputs["spring_l0"])
    l0_max = float(inputs["l0_max"])
    n_links_init = int(round(float(inputs["n_links_init"])))
    nlink_max = int(round(float(inputs["nlink_max"])))

    # Step through one exact floating-point split cycle matching src/filament.cpp:675-706
    lb = spring_l0
    steps = 0
    while True:
        steps += 1
        if lb + dL < l0_max:
            lb = lb + dL
        else:
            overshoot = (lb + dL) - l0_max
            break

    # In filament::update_length() (line 733): stops once get_nsprings() + 1 > nsprings_max
    # i.e., grows from n_links_init links up to nlink_max links via (nlink_max - n_links_init) splits
    n_splits_needed = max(0, nlink_max - n_links_init)
    delta_L_per_split = l0_max - spring_l0
    k_grow_nominal_um_s = kgrow * dL
    k_grow_realized_um_s = kgrow * (delta_L_per_split / float(steps)) if steps > 0 else 0.0
    mean_time = (n_splits_needed * steps) / kgrow if kgrow > 0 else float("inf")

    return {
        "k_grow_nominal": k_grow_nominal_um_s * 1000.0,
        "k_grow_realized": k_grow_realized_um_s * 1000.0,
        "steps_per_split": float(steps),
        "overshoot_per_split": max(0.0, overshoot),
        "mean_time_to_max_length": mean_time,
    }
