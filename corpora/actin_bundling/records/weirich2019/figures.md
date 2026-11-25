# Figure Inventory and Digitization Log — `weirich2019`

Complete page-by-page, panel-by-panel audit of every figure in `corpus/actin_bundling/2019_Weirich.pdf` (25 PDF pages: pp. 1–6 PNAS 116(23), 11125–11130 (2019); pp. 7–25 Supporting Information pp. SI-1–SI-19). All 38 quantitative figure readings across 16 quantitative panels were implemented as calibrated `mode: code` scripts (900–2400 dpi) under `figures/`, verified through `boxgraph figure check`, visually inspected via their magenta overlay PNGs in `figwork/`, and saved to `readings.json` (verified in bare-checkout mode with 38/38 readings loaded and 0 missing).

## Page 1 (printed p. 11125)
- No figures (Title, Abstract, Significance, and Introduction text only).

## Page 2 (printed p. 11126) — Figure 1: Molecular motors self-organize and deform biopolymer droplets
- **Fig. 1A**: Left out — experimental cartoon and fluorescence micrographs of filamin-cross-linked actin tactoids (`10 um`, `2 um` scale bars); qualitative morphology images with no quantitative axes.
- **Fig. 1B**: Left out — experimental cartoon of actin droplets prepared with myosin II filaments; no plot axes.
- **Fig. 1C (Upper & Lower)**: Left out — fluorescence micrographs of myosin II (white) centered at the midplane of actin tactoids (red) (`10 um`, `2 um` scale bars).
- **Fig. 1D (Upper & Lower)**: Left out — fluorescence micrographs of bisected and separated actin–myosin droplets (`10 um`, `2 um` scale bars).
- **Fig. 1E**: Left out — fluorescence time-lapse image sequence (`0', 4', 8', 12', 16', 20'`, `5 um` scale bar) of myosin puncta migrating, clustering, and bisecting a droplet (cluster dimensions are plotted and digitized in SI Fig. S1).
- **Fig. 1F**: Left out — fluorescence time-lapse image sequence (`0 s` to `80 s`, `2 um` scale bar) of a myosin cluster severing from actin droplets.
- **Fig. 1G (Upper — High % active motors morphology probability vs myosin concentration)**: Digitized all 3 marker series at 1200 dpi (`cal_res = (0.45, 0.43) px`, calibrating on ticks `[0, 2, 4, 6, 8, 10, 12] nM` and `[0, 0.2, 0.4, 0.6, 0.8, 1.0]`):
  - `fig1g_high_centered` (`d_fig1g_high_centered`, 3 red filled circles): `(2.46 nM, 0.8133), (6.17 nM, 0.1492), (12.26 nM, 0.0053)`.
  - `fig1g_high_bisected` (`d_fig1g_high_bisected`, 3 blue `x` markers): `(2.44 nM, 0.0009), (6.11 nM, 0.6574), (12.24 nM, 0.5907)`.
  - `fig1g_high_separated` (`d_fig1g_high_separated`, 3 purple open squares): `(2.44 nM, 0.0012), (6.12 nM, 0.1801), (12.24 nM, 0.3052)`.
  - Concentrations match Main p. 2 (`2.4 nM`), Main p. 3 (`6 nM` threshold), and Movie S2 (`12.3 nM`); probabilities sum to `0.815, 0.987, 0.901` (leaving 1–18% in uncategorized / myosin-free morphologies).
- **Fig. 1G (Lower — Lower % active motors morphology probability vs myosin concentration)**: Digitized all 3 marker series at 1200 dpi (`cal_res = (0.46, 0.47) px`):
  - `fig1g_low_centered` (`d_fig1g_low_centered`, 3 red filled circles): `(0.74 nM, 0.8086), (7.36 nM, 0.8230), (9.81 nM, 0.8658)`.
  - `fig1g_low_bisected` (`d_fig1g_low_bisected`, 3 blue `x` markers): `(0.72 nM, 0.0290), (7.34 nM, 0.0010), (9.79 nM, 0.0010)`.
  - `fig1g_low_separated` (`d_fig1g_low_separated`, 3 purple open squares): `(0.73 nM, 0.0010), (7.33 nM, 0.0010), (9.80 nM, 0.0010)`.
  - Concentrations match Movie S3/S6 (`0.7 nM`) and Movie S1 (`7.4 nM`); probabilities sum to `0.839, 0.825, 0.868`.

## Page 3 (printed p. 11127) — Figure 2: Motors migrate to regions of minimum nematic splay
- **Fig. 2A**: Left out — fluorescence time-lapse micrographs (`0', 2', 6', 12'`, `2 um` scale bar) of myosin centering after droplet coalescence.
- **Fig. 2B (Probability distribution of myosin II puncta vs normalized midplane distance `d_i / L_D`, `n = 50` droplets)**: Digitized both series at 1200 dpi (`cal_res = (0.45, 0.64) px`; notes a `+0.0023` (`+5.2 px`) rightward shift of the plotted data layer relative to the axis ticks):
  - `fig2b_t_lt_2min` (`d_fig2b_t_lt_2min`, 6 blue `x` markers): `0.3134, 0.3140, 0.2640, 0.0753, 0.0106, 0.0072` (sum `0.9845`).
  - `fig2b_t_gt_5min` (`d_fig2b_t_gt_5min`, 6 red filled circles): `0.6635, 0.2887, 0.0217, 0.0116, 0.0008, -0.0019` (sum `0.9844`).
  - **Finding (settles `q_fig2b_central_40_vs_60_percent`, logged in `q_digitized_fig2b_fig2e_histc_bin_shift`)**: At `t < 2 min`, the first two bins `[0, 0.4]` (central 40% of droplet) sum to `62.7%`, whereas the first three bins `[0, 0.6]` (central 60%) sum to `89.1% (~90%)`; at `t > 5 min`, the midplane bin holds `66.4%` (vs printed `65%`).
- **Fig. 2C (Main panel — multiple myosin puncta trajectories `d / 2L_D` vs time)**: Left out — ~6 densely overplotted trajectories of identical marker styles/colors that cross and tangle near `d / 2L_D = 0`, preventing unambiguous per-track separation (the representative single trajectory and fit are digitized below in the Inset).
- **Fig. 2C (Inset — single myosin II punctum trajectory `d (um)` vs time `0–1000 s` and exponential fit)**: Digitized at 2400 dpi (`cal_res = (0.69, 1.00) px`, inverted y-axis running from `-3` at top to `1` at bottom):
  - `fig2c_inset_data` (`d_fig2c_inset_data`): 94 red filled circles across `t = 11.5..960 s`, `d = -2.147..+0.419 um`.
  - `fig2c_inset_fit` (`d_fig2c_inset_fit`): 12 points along the solid red exponential curve (`t = 20..950 s`, `d = -1.522..+0.126 um`).
  - **Finding (settles `q_fig2c_fit_type_axis_and_unstated_params`, logged in `q_digitized_fig2c_inset_and_fig2d_speed_underprediction`)**: A 3-parameter single exponential `d(t) = -1.859 * exp(-t / 177.4 s) + 0.138 um` fits the plotted red curve across 1006 columns to RMS `0.0043 um`, giving `tau = 177.4 s = 2.96 min` (**3.4–6.8x faster** than Main p. 4 / SI p. 10's `10–20 min`) and initial speed `|d'(0)| = 10.5 nm/s`.
- **Fig. 2D (Initial myosin II puncta centering speed vs distance from midplane `d (um)`)**: Digitized all 3 series with full vertical error bars at 1200 dpi (`cal_res = (0.81, 0.34) px`, calibrating y strictly from ticks `[-2, 0, 2, 4, 6, 8, 10]` since the bottom box spine sits at `-2.20`, not `-2`; `error_bars: unstated`):
  - `fig2d_LD3p6` (`d_fig2d_LD3p6`, 2 black `x` markers): `(0.549, 0.187 [-0.329, 0.692]), (1.167, 0.740 [-1.822, 3.443]) nm/s`.
  - `fig2d_LD4p4` (`d_fig2d_LD4p4`, 8 blue open circles): `(0.902, 1.293 [0.958, 1.676]), (1.103, 1.396 [1.074, 1.755]), (1.531, 0.473 [0.309, 0.636]), (1.742, 1.293 [1.042, 1.559]), (1.792, 1.032 [0.944, 1.121]), (1.902, 5.294 [4.903, 5.723]), (2.182, 4.791 [4.503, 5.131]), (2.683, 8.792 [8.283, 9.342]) nm/s`.
  - `fig2d_LD5p1` (`d_fig2d_LD5p1`, 4 green filled diamonds): `(1.327, 1.713 [0.617, 2.814]), (1.559, 0.790 [-0.605, 2.198]), (1.608, 0.118 [-0.385, 0.836]), (2.769, 4.511 [3.844, 5.187]) nm/s`.
  - **Model validations (`rod_centering_force_velocity`, `F_0 = 32.5 pN*nm`, `Gamma = 15 pN*s/um`, `L_M = 1.5 um`)**:
    - `cmp_fig2d_v_LD3p6`: **`reproduced`** (2/2 `reproduced`, `max_abs = 0.245 nm/s`).
    - `cmp_fig2d_v_LD4p4`: **`not reproduced`** (1 `reproduced`, 4 `close`, 3 `not reproduced`, `max_abs = 7.36 nm/s`): predicts `0.65, 0.82, 1.42 nm/s` vs measured `5.29, 4.79, 8.79 nm/s` at `d = 1.90, 2.18, 2.68 um` (**6–8x too slow**).
    - `cmp_fig2d_v_LD5p1`: **`not reproduced`** (2 `reproduced`, 1 `close`, 1 `not reproduced`, `max_abs = 3.79 nm/s`): predicts `0.72 nm/s` vs measured `4.51 nm/s` at `d = 2.77 um` (**6.3x too slow**). Logged in `q_digitized_fig2c_inset_and_fig2d_speed_underprediction`.
- **Fig. 2E (Main histogram — alignment angle `theta_MD`, `n = 29` myosin)**: Digitized at 1200 dpi (`cal_res = (0.69, 0.79) px`):
  - `fig2e_alignment` (`d_fig2e_alignment`, 4 nonzero black bars): `(0.50 deg, 0.4830), (5.03 deg, 0.2756), (10.0 deg, 0.1723), (20.0 deg, 0.0690)` — exact integer counts **14, 8, 5, 2** out of 29 (sum `0.9999`).
  - **Finding (logged in `q_digitized_fig2b_fig2e_histc_bin_shift`)**: Exhibits the same `histc` left-bin-edge indexing error as Fig. 2B: Main p. 3's "> 90% within 10 deg" sums the first three bars (`27/29 = 93.1%`), which span `[0, 15)` deg for 5-deg left-edged bins, whereas `[0, 10)` deg (first two bars) contains only `22/29 = 75.9%`.
- **Fig. 2E (Inset)**: Left out — confocal fluorescence micrograph and z-projection of myosin embedded in a tactoid (`2 um`, `3 um`, `4 um` scale bars).
- **Fig. 2F (Continuum-model centering free energy `Delta F / k_B T` vs `d / L_D`)**: Digitized both solid curves at 1200 dpi (`cal_res = (0.40, 0.30) px`); embedded director-field cartoon left out:
  - `fig2f_LM1p5` (`d_fig2f_LM1p5`, 12 points along blue curve, `d/L_D = -0.6..+0.6`, `Delta F/k_B T = 0.231..8.046`).
  - `fig2f_LM0p33` (`d_fig2f_LM0p33`, 10 points along red curve, `d/L_D = -0.85..+0.85`, `Delta F/k_B T = 0.239..5.785`).
  - **Model comparisons (`rod_centering_energy_eq4_kBT`, `F_0 = 10 k_B T`)**:
    - With caption `L_T = 3.2 um`: `cmp_fig2f_LM1p5_caption` -> **`not reproduced`** (`max_abs = 1.31 k_B T`, +13–16% high); `cmp_fig2f_LM0p33_caption` -> **`not reproduced`** (`max_abs = 0.78 k_B T`, +12–16% high); `cmp_fig2f_LM0p30_caption` (`L_M = 0.3 um` from Main p. 4) -> `reproduced` within `0.3 ± 0.05` rounding (`max_abs = 0.135 k_B T`).
    - With experimental `L_D = 3.6 um` (Fig. 2D black `x`'s): `cmp_fig2f_LM1p5_LD3p6` -> **`reproduced`** (12/12, `max_abs = 0.017 k_B T`); `cmp_fig2f_LM0p33_LD3p6` -> **`reproduced`** (10/10, `max_abs = 0.010 k_B T`). Proves Fig. 2F was plotted with `L_D = 3.6 um`, not caption `L_T = 3.2 um` (`q_digitized_LD3p6_figs_2f_S4_S5`).

## Page 4 (printed p. 11128) — Figure 3: Motor clusters locally distort nematic alignment
- **Fig. 3A**: Left out — fluorescence micrographs of actin tactoids with myosin puncta and clusters (`2 um` scale bars); qualitative imaging.
- **Fig. 3B**: Left out — fluorescence micrographs of thin crowded actin nematic layers (`< 1 um`) with myosin clusters (`5 um`, `2 um` scale bars); qualitative imaging.

## Page 5 (printed p. 11129) — Figure 4: Molecular motors bisect biopolymer droplets
- **Fig. 4A**: Left out — theoretical two-droplet free energy vs `L_A / (L_A + L_B)` for `w/gamma = 1.5` (red) and `3.5` (blue), whose y-axis ("Free Energy") has 7 unnumbered ticks and no numerical scale ("normalized to display minima"). The numbered-y-axis counterpart of these Eq. 7 curves is digitized in SI Fig. S9.
- **Fig. 4B (Measured daughter-droplet lengths `L_B` vs `L_A` in bisected droplets)**: Digitized at 1200 dpi (`cal_res = (0.38, 0.36) px`):
  - `fig4b_bisected_lengths` (`d_fig4b_bisected_lengths`): all 46 black filled circles across `L_A = 3.06..7.85 um`, `L_B = 3.08..7.94 um`. OLS regression gives `L_B = 1.0038 * L_A - 0.0390 um` (`R^2 = 0.9606`), with `RMS(L_B - L_A) = 0.254 um`, `mean |L_B - L_A| = 0.198 um`, `max |L_B - L_A| = 0.713 um`.
  - **Model validation (`equal_droplet_division_length`, `L_B = L_A`)**: `cmp_fig4b_equal_length` -> **`not reproduced`** point-by-point (2 `reproduced`, 4 `close`, 40 `not reproduced`, `max_abs = 0.713 um`), reflecting biological droplet-to-droplet scatter (~0.25 um RMS) around the exact slope-1 trend relative to the `±0.02 um` marker reading precision.
- **Fig. 4C**: Left out — fluorescence time-lapse micrographs (`-0.25'` to `9'`, `2 um` scale bar) of an isolated droplet coalescing with a bisected droplet.
- **Fig. 4D (Phase space `V_droplet / V_cluster` vs `w/gamma`, Eq. 7)**: Digitized both dot-sequence boundaries at 1200 dpi (`cal_res = (0.44, 0.05) px`, `log10` y-axis calibrated on ticks `[1, 2, 5, 10, 20]`):
  - `fig4d_equal_boundary` (`d_fig4d_equal_boundary`, 9 points along black dot sequence, `w/gamma = 2.1..5.8`, `V_d/V_c = 1.517..13.35`).
  - `fig4d_unequal_boundary` (`d_fig4d_unequal_boundary`, 7 points along red dot sequence, `w/gamma = 2.6..3.8`, `V_d/V_c = 3.334..23.79`).
  - **Model comparisons**: `cmp_fig4d_equal` (`spherical_droplet_equal_boundary`) -> **`reproduced`** (9/9, `max_abs = 0.010`); `cmp_fig4d_unequal` (`spherical_droplet_unequal_boundary`) -> **`reproduced`** (7/7, `max_abs = 0.054`). Both numerical root-finder refinements (`ref_fig4d_equal_tol`, `ref_fig4d_unequal_tol`) settle to `< 2.2e-11` (`belongs_to: model`).
- **Fig. 4E**: Left out — cartoon schematic summarizing passive vs active transitions in actomyosin droplets (panels `i–v`, with panel iv mislabeled `vi`).

## Page 6 (printed p. 11130)
- No figures (Materials and Methods, Acknowledgments, and References text only).

## Pages 7–12 (printed pp. SI-1–SI-6)
- No figures (SI Appendix cover page, Experimental Methods, and Theory Sections 1–3 derivations).

## Page 13 (printed p. SI-7) — SI Figure S1: Myosin cluster width vs length during clustering
- **Fig. S1 (Myosin cluster width vs length in Fig. 1E / Movie S3)**: Digitized all 3 time-point marker classes at 900 dpi (`cal_res = (0.54, 0.95) px`); qualitative shaded "puncta" vs "cluster" dashed boundary left out:
  - `figS1_0min` (`d_figS1_0min`, 6 black `x` markers): `(0.731, 0.336), (0.798, 0.336), (0.998, 0.536), (0.998, 0.602), (1.464, 0.736), (1.598, 0.803) um`.
  - `figS1_8min` (`d_figS1_8min`, 3 blue filled circles): `(1.066, 0.599), (1.600, 0.999), (1.867, 0.932) um`.
  - `figS1_20min` (`d_figS1_20min`, 1 red open diamond): `(1.800, 1.202) um`.

## Page 14 (printed p. SI-8) — SI Figure S2
- **Fig. S2**: Left out — fluorescence micrographs of low-%-active myosin at `17 nM` and `58 nM` in a thin actin nematic layer (scale bar misprinted `10 uM` instead of `10 um`); no quantitative axes.

## Page 15 (printed p. SI-9) — SI Figure S3
- **Fig. S3A & S3B**: Left out — geometric schematics of a rod embedded in a bipolar tactoid and the intersecting circular-arc director construction; no quantitative axes.

## Page 16 (printed p. SI-10) — SI Figure S4: Dimensionless rod centering free energy vs `x_c / L_D`
- **Fig. S4 (Shifted Eq. 4 free energy `Delta F` vs `x_c / L_D`)**: Digitized all 3 solid curves at 900 dpi (`cal_res = (0.40, 0.30) px`):
  - `figS4_ratio0p4` (`d_figS4_ratio0p4`, blue, 10 pts `x_c/L_D = -0.7..+0.7`, `Delta F = 0.094..1.380`).
  - `figS4_ratio0p2` (`d_figS4_ratio0p2`, red, 8 pts `x_c/L_D = -0.8..+0.8`, `Delta F = 0.065..0.985`).
  - `figS4_ratio0p1` (`d_figS4_ratio0p1`, green, 8 pts `x_c/L_D = -0.9..+0.9`, `Delta F = 0.045..0.943`).
  - **Model comparisons (settles `q_figS4_figS5_caption_vs_legend_ratios`, logged in `q_digitized_LD3p6_figs_2f_S4_S5`)**:
    - Caption `L_D = 3 um` (`cmp_figS4_{0p4,0p2,0p1}_caption`): all 3 **`not reproduced`** (`max_abs = 0.594, 0.263, 0.231`, +22–43% high).
    - Legend rounded ratios `0.4, 0.2, 0.1` (`cmp_figS4_{0p4,0p2,0p1}_legend`): `reproduced` only within 1-digit `±0.05` rounding (`max_abs = 0.088, 0.018, 0.099`, systematic 4–11% error on blue and green).
    - Actual `L_D = 3.6 um` (`L_M/L_D = 0.4167, 0.2028, 0.0917`; `cmp_figS4_{0p4,0p2,0p1}_LD3p6`): all 3 **`reproduced`** (26/26 points, `max_abs = 0.0018, 0.0013, 0.0050`, 20–50x smaller error).

## Page 17 (printed p. SI-11) — SI Figure S5: Transverse rod centering free energy vs `y_c / L_y`
- **Fig. S5 (Exact Eq. 3 LHS transverse integral `Delta F / f_0` in `um` vs `y_c / L_y` at `x_c = 1 um`)**: Digitized all 3 solid curves at 900 dpi (`cal_res = (0.00, 0.39) px`):
  - `figS5_ratio0p4` (`d_figS5_ratio0p4`, blue, 10 pts `y_c/L_y = -0.9..+0.9`, `Delta F/f_0 = 0.0317..0.2526 um`).
  - `figS5_ratio0p2` (`d_figS5_ratio0p2`, red, 10 pts `y_c/L_y = -0.9..+0.9`, `Delta F/f_0 = 0.0148..0.1195 um`).
  - `figS5_ratio0p1` (`d_figS5_ratio0p1`, green, 10 pts `y_c/L_y = -0.9..+0.9`, `Delta F/f_0 = 0.0066..0.0538 um`).
  - **Model comparisons (`rod_centering_energy_transverse`)**:
    - Caption `L_D = 3 um, L_y = 1 um` (`cmp_figS5_{0p4,0p2,0p1}_caption`): all 3 **`not reproduced`** (`max_abs = 0.065, 0.028, 0.013 um`, +26–30% high).
    - Consistent `L_D = 3.6 um, L_y = 1.2 um` (aspect ratio 3; `cmp_figS5_{0p4,0p2,0p1}_LD3p6`): all 3 **`reproduced`** (30/30 points, `max_abs = 0.00012, 0.00007, 0.00008 um`). Quadrature refinement `ref_figS5_transverse_nquad` settles to `8.9e-16` (`belongs_to: model`).

## Page 18 (printed p. SI-12) — SI Figure S6: Harmonic force constant `k / F_0` vs `L_D`
- **Fig. S6 (Eq. 5 harmonic spring constant `k / F_0` vs `L_D`)**: Digitized both solid curves at 900 dpi (`cal_res = (0.40, 0.39) px`, calibrating y strictly from the 7 right-spine ticks `[0.2..1.4]` since the bottom spine sits at `0.112`, not `0`):
  - `figS6_LM1p5` (`d_figS6_LM1p5`, blue, 9 pts `L_D = 2.1..5.9 um`, `k/F_0 = 1.346..0.1253 um^-2`).
  - `figS6_LM0p33` (`d_figS6_LM0p33`, red — caption misprints "green", 9 pts `L_D = 2.1..5.9 um`, `k/F_0 = 0.9251..0.1179 um^-2`).
  - **Model comparisons (`rod_centering_spring_constant`, `k_over_F0`)**: `cmp_figS6_LM1p5` -> **`close`** (8 `reproduced`, 1 `close` at `L_D = 2.1 um`, `max_abs = 0.0047 um^-2` = 0.35%); `cmp_figS6_LM0p33` -> **`reproduced`** (9/9, `max_abs = 0.0027 um^-2`).

## Page 19 (printed p. SI-13) — SI Figure S7: Overdamped centering velocity `V` vs `x_c / L_D`
- **Fig. S7**: Left out — theoretical velocity curves for `L_D = 4, 5, 6` whose y-axis (`V`) carries zero tick marks, numbers, origin, or units, making quantitative y-calibration impossible. Independently calibrating the numbered x-axis (`0.2..0.8`, `cal_res_x = 0.00 px`) locates the top-frame asymptote crossings at `x_c / L_D = 0.8080, 0.8467, 0.8724`, quantitatively confirming the unstated `L_M = 1.5 um` (`1 - L_M/(2 L_D) = 0.8125, 0.8500, 0.8750`).

## Page 20 (printed p. SI-14) — SI Figure S8
- **Fig. S8A & S8B**: Left out — geometric schematics of one and two bipolar tactoids adhered orthogonally to a spherical colloid of radius `a`; no plot axes.

## Page 21 (printed p. SI-15) — SI Figure S9: Energetics of spheroidal tactoid division (Eq. 7)
- **Fig. S9 (Free energy `F` vs `V_1 / (V_1 + V_2)`, Eq. 7 at unstated `V_d/V_c = 4`)**: Digitized all 3 curves at 900 dpi (`cal_res = (2.15, 1.50) px` on the embedded 300-dpi Mathematica bitmap, i.e., `< 0.7` native bitmap pixels; calibrating x strictly on the 21 padded ticks `[0.0..1.0]`):
  - `figS9_w1p5` (`d_figS9_w1p5`, solid purple, 11 pts `0.05..0.95`, `F = 2.179..2.340`).
  - `figS9_w2p5` (`d_figS9_w2p5`, dashed red, 13 pts `0.05..0.95`, asymmetric minima `2.071` at `0.15` & `0.85`, central maximum `2.133` at `0.5`).
  - `figS9_w3p5` (`d_figS9_w3p5`, solid black, 11 pts `0.05..0.95`, central minimum `2.026` at `0.5`).
  - **Model comparisons (`spherical_droplet_division_energy`, raw Eq. 7 at `V_d/V_c = 4`)**: `cmp_figS9_w1p5`, `cmp_figS9_w2p5`, `cmp_figS9_w3p5` -> all 3 **`not reproduced`** (35/35 points, `max_abs = 0.300, 0.852, 1.500`), revealing the exact undocumented affine transformations `F_S9 = F_Eq7 + 0.3` (`w/gamma = 1.5`), `4.0*F_Eq7 - 3.0` (`w/gamma = 2.5`), and `F_Eq7 + 1.5` (`w/gamma = 3.5`) (`q_digitized_figS9_S10_exact_affine_rescaling`).

## Page 22 (printed p. SI-16) — SI Figure S10: Energetics of highly elongated tactoid division (Eq. 8)
- **Fig. S10 (Free energy `F` vs `V_1 / (V_1 + V_2)`, Eq. 8 at unstated `V_d/V_c = 1`)**: Digitized all 3 curves at 900 dpi (`cal_res = (0.25, 0.23) px`, calibrating x strictly on the 21 padded ticks `[0.0..1.0]`):
  - `figS10_w0p1` (`d_figS10_w0p1`, solid blue, 11 pts `0.05..0.95`, `F = 1.763..2.091`).
  - `figS10_w0p5` (`d_figS10_w0p5`, dashed red, 11 pts `0.05..0.95`, `F = 0.697..0.947`).
  - `figS10_w2p0` (`d_figS10_w2p0`, solid black, 11 pts `0.05..0.95`, central minimum `0.4885` at `0.5`).
  - **Model comparisons (`elongated_tactoid_division_energy`, raw Eq. 8 at `V_d/V_c = 1`)**: `cmp_figS10_w0p1`, `cmp_figS10_w0p5`, `cmp_figS10_w2p0` -> all 3 **`not reproduced`** (33/33 points, `max_abs = 0.945, 0.498, 2.718`), revealing the exact undocumented affine transformations `F_S10 = 2.0*F_Eq8 - 0.2` (`W_el = 0.1`), `3.0*F_Eq8 - 0.4` (`W_el = 0.5`), and `1.3*F_Eq8 + 3.3` (`W_el = 2.0`), matching to RMS `<= 0.00013` (`q_digitized_figS9_S10_exact_affine_rescaling`).

## Page 23 (printed p. SI-17) — SI Figure S11: Phase map of elongated tactoid division (Eq. 8)
- **Fig. S11 (Phase map `V_d / V_c` vs `w a^2 gamma / K^2`, Eq. 8)**: Digitized both dot-sequence boundaries at 900 dpi (`cal_res = (0.30, 0.28) px`, `log10` y-axis calibrated on ticks `[1, 2, 5, 10, 20]`):
  - `figS11_equal_boundary` (`d_figS11_equal_boundary`, 9 points along black dot sequence, `W_el = 1.3..4.8`, `V_d/V_c = 1.399..36.65`).
  - `figS11_unequal_boundary` (`d_figS11_unequal_boundary`, 8 points along red dot sequence, `W_el = 1.75..2.40`, `V_d/V_c = 3.327..32.69`).
  - **Model comparisons**: `cmp_figS11_equal` (`elongated_tactoid_equal_boundary`, `2*(2*W_el/3)^(5/2)`) -> **`reproduced`** (9/9, `max_abs = 0.010`); `cmp_figS11_unequal` (`elongated_tactoid_unequal_boundary`, `1 + (W_el/(3 - W_el))^(5/2)`) -> **`reproduced`** (8/8, `max_abs = 0.310` = 0.95% at `V_d/V_c = 32.7`), confirming the implicit `beta = 1` assumption noted in `q_eq8_missing_length_ratio_and_large_alpha_figS11`.

## Pages 24–25 (printed pp. SI-18–SI-19) — Supporting Movies S1–S8 and References
- **Movies S1–S8 captions (p. 24)** and **SI References (p. 25)**: No figure panels (text only).
