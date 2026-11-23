# Figure Readings — `claessens2008`

Calibrated `mode: code` digitization of every quantitative plot across all 7 PDF pages of Claessens, Semmrich, Ramos & Bausch, *Proc. Natl. Acad. Sci. USA* **105**(26), 8819–8822 (2008) plus Supporting Information (`corpus/actin_bundling/2008_Claessens.pdf`, pp. 1–4 main article, pp. 5–7 SI). Main-text panels were rendered at **1200 dpi** (plot boxes 1041–2849 px) and SI Fig. S1 at **900 dpi** (plot box 4266 × 3462 px). All axes were calibrated directly from detected frame borders and inward ticks via `figures/_calib.py`, and all **22 overlay PNGs** in `figwork/` were generated with `boxgraph figure check` and visually inspected to confirm every magenta cross (and vertical error-bar span) sits cleanly on its target marker, curve, or error bar.

- **Page 1 (8819) — No figures:** Title, abstract, and main text only.

- **Fig. 1A (p. 2 / 8820) — Left out:**
  - *Reason:* Fluorescence micrograph of TRITC phalloidin-labeled F-actin/fascin bundles ($c_a = 0.04\text{ mg/ml}$, $R^* = 1$, scale bar 10 µm); micrograph with no plotted axes or quantitative curve.

- **Fig. 1B (p. 2 / 8820) — Left out:**
  - *Reason:* Fluorescence micrograph at $c_a = 0.1\text{ mg/ml}$, $R^* = 1$ (scale bar 10 µm); micrograph with no plotted axes or quantitative curve.

- **Fig. 1C (p. 2 / 8820) — Read (`fig1c_diameter_data`, `fig1c_powerlaw_fit`); Inset left out:**
  - *Panels & series:* Log-log plot of TEM bundle diameter $D$ (nm) vs. effective bound ratio $R^*$ at $c_a = 0.1\text{ mg/ml}$:
    1. `fig1c_diameter_data` (`D_fig1c_measured`, `origin: measured`): all 8 filled black circular markers with unstated vertical error bars (`error_bars: unstated`; caps traced directly from vertical stems and horizontal caps):
       - $R^* = 0.00890$: $D = 8.41\text{ nm}$ (error bar `[7.72, 9.14]`)
       - $R^* = 0.01778$: $D = 11.05\text{ nm}$ (`[10.42, 11.74]`)
       - $R^* = 0.04434$: $D = 14.39\text{ nm}$ (`[13.61, 15.21]`)
       - $R^* = 0.08824$: $D = 16.65\text{ nm}$ (`[15.03, 18.27]`)
       - $R^* = 0.1739$: $D = 19.94\text{ nm}$ (`[18.36, 21.54]`)
       - $R^* = 0.4125$: $D = 27.12\text{ nm}$ (`[25.35, 28.87]`)
       - $R^* = 0.7023$: $D = 28.14\text{ nm}$ (`[27.15, 29.22]`)
       - $R^* = 0.9003$: $D = 28.04\text{ nm}$ (`[26.92, 29.22]`)
    2. `fig1c_powerlaw_fit` (`D_fig1c_powerlaw_fit`, `origin: fitted`): dashed straight power-law line sampled at 8 interior $R^*$ points (`0.008..0.8`), fitted from 20 clean isolated dash centroids (`max_res = 0.56 px`, log-log slope `0.3269`, $D(R^*=1) = 38.00\text{ nm}$).
    3. *Fig. 1C Inset — Left out:* TEM micrograph of F-actin/fascin bundles (no plotted axes).
  - *Calibration (1200 dpi):* Log10 x-axis from 2 interior major ticks + right border (`0.01, 0.1, 1`, residual `0.83 px`, left border `0.00500`, all 22 minor ticks verified); log10 y-axis from 6 matched left/right ticks (`7, 8, 9, 10, 20, 30`, residual `0.87 px`, frame `[6.00, 40.02]`).
  - *Verdicts & findings (`q_fig1c_figS1_slope_and_Dmax_anchor`):*
    - `cmp_D_scaling_fig1c_fit` (`reproduction`, `diameter_scaling_Rstar`, fitted prefactor `35.56 nm`): **not reproduced** (2 close, 6 not reproduced) because the drawn dashed line has slope `0.3269` rather than the text's 1-sig-fig rounded `0.3` (`8.35` vs. `7.84 nm` at $R^* = 0.008$; `33.25` vs. `35.33 nm` at $R^* = 0.8$).
    - `cmp_D_scaling_fig1c_data` (`validation`, `diameter_scaling_Rstar`, fitted prefactor `33.92 nm`): **not reproduced** (4 reproduced, 3 close, 1 not reproduced at $R^* = 0.90$): captures the rising regime through $R^* \approx 0.41$ while departing at the saturated plateau $D \approx 28.1\text{ nm}$ ($R^* = 0.70, 0.90$), confirming `q_tem_vs_saxs_bundle_diameter`.

- **Fig. 2A, 2B, 2C (p. 2 / 8820) — Left out:**
  - *Reason:* Fluorescence micrographs of confined F-actin/fascin bundles in emulsion droplets at $R^* = 1$ (scale bars 10 µm), carrying bottom-left color-key circles (A red = single ring, B green = two bundles/split ring, C dark blue = complicated structures); no plotted axes.

- **Fig. 2D (p. 2 / 8820) — Left out:**
  - *Reason:* TEM micrograph of an actin bundle extracted from confined rings (scale bar 20 nm, ~5–6 filaments across); micrograph with no plotted axes.

- **Fig. 2E (p. 2 / 8820) — Read (`fig2e_nf20_boundary`, `fig2e_single_ring_max_Dd`, `fig2e_two_bundle_min_Dd`):**
  - *Series:* Log-log confinement phase diagram $c_a$ (mg/ml) vs. droplet diameter $D_d$ (µm) at $R^* = 1$:
    1. `fig2e_nf20_boundary` (`ca_fig2e_nf20_boundary`, `origin: predicted`): thick dashed black line ("$n_f = 20$") fitted across 15 clean dashes and sampled at 8 droplet diameters (`6.5..23 µm`). Measures log-log slope `-2.014` and uniform product $c_a D_d^2 = 6.43 \pm 0.04\,(\text{mg/ml})\cdot\mu\text{m}^2$.
    2. `fig2e_single_ring_max_Dd` (`ca_fig2e_single_ring_max`, `origin: measured`): rightmost red single-ring circle center at each of the 5 actin concentrations (`c_a = 0.10, 0.05, 0.04, 0.03, 0.02 mg/ml` $\to D_d = 6.88, 10.74, 16.08, 17.50, 22.71\,\mu\text{m}$; $c_a D_d^2 = 4.73\text{–}10.25$).
    3. `fig2e_two_bundle_min_Dd` (`ca_fig2e_two_bundle_min`, `origin: measured`): leftmost green two-bundle circle center at each concentration (`0.10, 0.05, 0.04, 0.03, 0.02 mg/ml` $\to D_d = 9.60, 13.12, 15.73, 26.96, 25.96\,\mu\text{m}$; $c_a D_d^2 = 8.53\text{–}21.74$).
    - *Note on overlapping interior circles:* Individual interior droplets in the dense horizontal clusters overlap heavily across six rainbow transition colors (red, orange, yellow, green, cyan, dark blue) and cannot be uniquely deconvolved; the upper single-ring and lower two-bundle edges cleanly bracket the experimental phase boundary at every measured concentration.
  - *Calibration (1200 dpi):* Log10 x-axis from 16 matched top/bottom ticks (`4..10, 20..100`, residual `2.35 px` across 2284 px = `0.10%`; frame `[3.01, 110.1] µm`); log10 y-axis from 9 interior ticks + bottom border (`0.01..0.10`, residual `2.32 px` across 1760 px = `0.13%`; top border `0.200 mg/ml`).
  - *Verdicts & settled open question (`q_droplet_ring_formula_and_fig2e_factor`, `q_fig2e_quantified_factor_2p22`):*
    - `cmp_droplet_crit_ca_fig2e_line` (`reproduction`, `droplet_critical_concentration` with $n_f = 20$, $M_{\text{actin}} = 42\text{ kDa}$, $d = 2.9\text{ nm}$): **not reproduced** (0/8): model predicts $c_a D_d^2 = 2.886\,(\text{mg/ml})\cdot\mu\text{m}^2$, uniformly **2.23× lower** than the plotted dashed line (`6.43`; still `2.11×` with $d = 2.75\text{ nm}$).
    - `cmp_droplet_crit_ca_fig2e_single_ring` (`validation`): **not reproduced** (0/5, model `1.64×–3.55×` below measured rightmost red single-ring droplets).

- **Fig. 3A & Inset (p. 3 / 8821) — Read (`fig3a_saxs_R0p5`, `fig3a_saxs_R0p05`, `fig3a_saxs_R0p005`, `fig3a_inset_R0p5`):**
  - *Series:* Circularly averaged SAXS intensity $I(q)$ (a.u.) vs. $q$ ($\text{nm}^{-1}$) at $c_a = 2\text{ mg/ml}$:
    1. `fig3a_saxs_R0p5` (`I_saxs_fig3a_R0p5`, `origin: measured`): red curve, 18 points on $q \in [0.35, 1.95]\text{ nm}^{-1}$, capturing the sharp `(1,0)` peak at `q = 0.5852 nm^-1` (matching printed `0.585 nm^-1` to `0.03%`) and peaks `a`, `b`.
    2. `fig3a_saxs_R0p05` (`I_saxs_fig3a_R0p05`, `origin: measured`): blue curve, 13 points on $q \in [0.35, 1.95]\text{ nm}^{-1}$.
    3. `fig3a_saxs_R0p005` (`I_saxs_fig3a_R0p005`, `origin: measured`): green curve, 12 points on $q \in [0.35, 1.95]\text{ nm}^{-1}$.
    4. `fig3a_inset_R0p5` (`I_saxs_fig3a_inset_R0p5`, `origin: measured`): magnified $R^* = 0.5$ red curve in the Inset, 18 points on $q \in [0.85, 1.92]\text{ nm}^{-1}$ spanning all labeled features.
  - *Calibration (1200 dpi):* Main linear x-axis (`0.4..2.0`, 5 marks, residual `2.30 px` / 2250 px) and log10 y-axis (`0.005..0.05`, 10 ticks, residual `1.79 px` / 2849 px); Inset linear x-axis (`0.8..2.0`, 7 marks, residual `1.44 px` / 1244 px) and **log10** y-axis (`0.50..0.63`, 14 marks, residual `0.73 px` vs. linear `20.44 px` — confirming logarithmic ordinate).
  - *Settled open question (`q_fig3a_inset_12_21_degeneracy`, `q_fig3a_fig3d_calibrated_verification`):* Calibrated Inset peak/arrow positions confirm `(1,1)` arrow `1.051` (vs. $\sqrt{3}q_{10} = 1.013$), `(2,0)/a` peak `1.198` (printed `1.20`), `b` peak `1.341` / arrow `1.344` (printed `1.35`), `c` peak `1.456` / arrow `1.454` (printed `1.46`), `(2,1)` bump `1.561` / arrow `1.585` (vs. $\sqrt{7}q_{10} = 1.548$), and `(3,0)` bump `1.743` / arrow `1.765` (vs. $3q_{10} = 1.755$), whereas the `(1,2)` arrow points to **`q = 1.693 nm^-1`** on a featureless descending slope with no peak and no 2D hexagonal Bragg reflection, decisively confirming the crystallographic mislabeling. (No runnable SAXS intensity model exists in the record.)

- **Fig. 3B (p. 3 / 8821) — Left out:**
  - *Reason:* 2D SAXS detector diffraction image ($R^* = 0.5$, $c_a = 2\text{ mg/ml}$) with schematic $q_z, q_r$ arrows only and no numerical axis ticks or color-bar scale.

- **Fig. 3C (p. 3 / 8821) — Read (`fig3c_fit_R0p05`, `fig3c_fit_R0p1`, `fig3c_fit_R0p5`, `fig3c_data_R0p05`, `fig3c_data_R0p1`, `fig3c_data_R0p5`):**
  - *Series:* Angularly averaged axial SAXS profiles $I(q_z)$ (a.u.) vs. $q_z$ ($\text{nm}^{-1}$):
    - Continuous four-sphere model fit curves (`origin: fitted`, isolated from hollow marker outlines via radius-4 disk morphological opening): `fig3c_fit_R0p05` (red, 20 points, peak `I = 0.246` at `1.18`), `fig3c_fit_R0p1` (blue, 18 points, peaks `0.454` at `1.18` and `0.118` at `1.32`), `fig3c_fit_R0p5` (green, 15 unobstructed points, peak `0.694` at `1.19` and secondary rise `0.101` at `1.32`; `1.34–1.36` omitted where occluded by crossing markers).
    - Measured marker profiles (`origin: measured`): `fig3c_data_R0p05` (red open circles, 17 points), `fig3c_data_R0p1` (blue open squares, 16 points), `fig3c_data_R0p5` (green triangles, 17 points — note caption calls them "filled triangle" though rendered hollow; clearly shows wings staying at `I ≈ 0.03–0.06` where the continuous fit curves drop to `0.00`).
  - *Calibration (1200 dpi):* Linear x-axis from 4 ticks (`1.1, 1.2, 1.3, 1.4`, residual `0.15 px`), linear y-axis from 7 matched left/right ticks (`0.1..0.7`, residual `0.36 px`, $I = 0$ at row `1084.4`).
  - *Model note:* `four_sphere_diffraction_model` is non-runnable (`q_four_sphere_model_unstated`), so no model comparison is run.

- **Fig. 3D (p. 3 / 8821) — Read (`fig3d_overtwist_ca0p8`, `fig3d_overtwist_ca2`, `fig3d_overtwist_ca3`):**
  - *Series (`origin: measured`):* All 12 measured helical overtwist points $\Delta\theta$ ($^\circ$) vs. $R^*$, quantized to `0.1 deg`:
    - `fig3d_overtwist_ca0p8` ($c_a = 0.8\text{ mg/ml}$, 4 green triangles, geometric centroids): `(0.04920, 0.601°)`, `(0.09832, 0.901°)`, `(0.4852, 1.001°)`, `(0.8824, 1.001°)`.
    - `fig3d_overtwist_ca2` ($c_a = 2\text{ mg/ml}$, 4 red squares, bounding-box centers recovering occluded interiors): `(0.02983, 0.400°)`, `(0.04979, 0.700°)`, `(0.09959, 0.900°)`, `(0.4957, 1.100°)`.
    - `fig3d_overtwist_ca3` ($c_a = 3\text{ mg/ml}$, 4 orange circles, centroids): `(0.02987, 0.300°)`, `(0.04966, 0.700°)`, `(0.09930, 1.000°)`, `(0.4939, 1.100°)`.
  - *Calibration (1200 dpi):* Log10 x-axis from 3 major ticks (`0.01, 0.1, 1`, residual `0.33 px`; all 19 log ticks have residual `0.41 px`); linear y-axis from 6 ticks (`0.0..1.0`, residual `0.33 px`).
  - *Settled / informed open questions (`q_Rstar_stoichiometry_contradiction`, `q_twist_energy_missing_half_and_constancy`, `q_twist_radian_unit_and_high_Rstar`, `q_fig3a_fig3d_calibrated_verification`):*
    - All 12 $R^*$ abscissas match the 1:1 quadratic `effective_bound_ratio` formula at $R \in \{0.03, 0.05, 0.1, 0.5, 1.0\}$ ($K_d = 0.3\,\mu\text{M}$, $M_{\text{actin}} = 42\text{ kDa}$) to `< 0.4%` (up to $R^* = 0.882$ at $c_a = 0.8\text{ mg/ml}, R = 1$), confirming `q_Rstar_stoichiometry_contradiction`.
    - Evaluating $E_{\text{twist}}/(R^* k_B T)$ at all 12 points with $\tau = 8\times 10^{-26}\text{ N}\cdot\text{m}^2$, $d = 2.9\text{ nm}$, $T = 20^\circ\text{C}$ directly tests constancy:
      - Printed formula (no `1/2`): `15.23, 17.14, 4.29, 2.36 k_B T` ($0.8\text{ mg/ml}$); `11.13, 20.42, 16.88, 5.07 k_B T` ($2\text{ mg/ml}$); `6.24, 20.49, 20.91, 5.08 k_B T` ($3\text{ mg/ml}$) — varies 9-fold and **exceeds $\Delta G \approx 15\,k_B T$ at all 6 points with $R^* = 0.05\text{–}0.10$**.
      - Linear-elastic `1/2` variant: `7.62, 8.57, 2.14, 1.18 k_B T` ($0.8\text{ mg/ml}$); `5.56, 10.21, 8.44, 2.53 k_B T` ($2\text{ mg/ml}$); `3.12, 10.25, 10.46, 2.54 k_B T` ($3\text{ mg/ml}$) — peaks right at `~10 k_B T` and stays strictly `< 15 k_B T` everywhere, though still non-constant across $R^*$.

- **Fig. 4A (p. 3 / 8821) — Left out:**
  - *Reason:* Fluorescence micrograph of composite actin/fascin/espin bundles (scale bar 5 µm); micrograph with no plotted axes.

- **Fig. 4B & Inset (p. 3 / 8821) — Left out:**
  - *Reason:* Main panel is a fluorescence micrograph of frayed actin/fascin/$\alpha$-actinin bundles. The overlaid Inset shows a single qualitative fluorescence line-scan trajectory along the frayed bundle with calibrated y-ticks (`300–500` a.u.) and 5 horizontal dashed guide lines (~262, 294, 358, 424, 490 a.u., corresponding to 1 + 2 + 2 + 2 ≈ 7 sub-bundle steps of ~32 a.u.), but has **no ticks, numbers, or units on its horizontal axis** ("position along the bundle axis"), so it cannot be calibrated in x.

- **Fig. 4C (p. 3 / 8821) — Left out:**
  - *Reason:* TEM micrograph of a WT human espin/actin bundle ($c_a = 0.95\,\mu\text{M}$, $R_{\text{espin}} = 4$, scale bar 50 nm); micrograph with no plotted axes.

- **Fig. 4D (p. 3 / 8821) — Left out:**
  - *Reason:* TEM micrograph of a thick fascin + espin/actin bundle ($c_a = 0.95\,\mu\text{M}$, $R_{\text{espin}} = 4$, $R_{\text{fascin}} = 1$, scale bar 50 nm); micrograph with no plotted axes.

- **Page 4 (8822) — No figures:** Conclusions, Methods, Acknowledgments, and References only.

- **Page 5 (SI 1) — No figures:** SI Methods text and bundle-geometry equation only.

- **Fig. S1 (p. 6 / SI 2) — Read (`figS1_blue_surface`, `figS1_green_inside`, `figS1_experiment`, `figS1_powerlaw`):**
  - *Series:* Log-log plot of cross-linker/actin ratio ($R^*$ ordinate) vs. bundle filament number $n$:
    1. `figS1_blue_surface` (`R_figS1_blue_surface`, `origin: predicted`, y = `molar_ratio_abp_actin`): blue solid theoretical curve (surface occupation included, $u = 13$), 11 points on $n \in [1.5, 80]$ (avoiding the red dashed crossing near $n \approx 17.5$).
    2. `figS1_green_inside` (`R_figS1_green_inside`, `origin: predicted`, y = `molar_ratio_abp_actin`): green solid theoretical curve (inside occupation only, $u = 13$), 12 points on $n \in [1.6, 80]$ (avoiding the red dashed crossing near $n \approx 10.6$).
    3. `figS1_experiment` (`Rstar_figS1_experiment`, `origin: measured`, y = `effective_bound_ratio`): all 8 red filled experimental circles (isolated from the black-outlined schematic hexagons): `(1.780, 0.00888)`, `(3.071, 0.01769)`, `(5.189, 0.04422)`, `(6.963, 0.08758)`, `(9.962, 0.1732)`, `(18.51, 0.4113)`, `(19.90, 0.7024)`, `(19.73, 0.8943)`.
    4. `figS1_powerlaw` (`Rstar_figS1_powerlaw`, `origin: fitted`, y = `effective_bound_ratio`): red dashed scaling line fitted across 33 clean dashes (`fit_max_res = 3.21 px`, log-log slope `1.5627`), sampled at 7 points (`n = 2..15`).
  - *Calibration (900 dpi):* Log10 x-axis from left/right borders + 17 interior ticks (`1..10, 20..100`, residual `2.75 px` / 4266 px = `0.06%`; 3 major ticks residual `1.33 px`); log10 y-axis from 23 ticks (`0.006..1.0`, residual `0.79 px` / 3462 px = `0.02%`; frame `[0.00500, 2.000]`).
  - *Verdicts & major findings (`q_figS1_surface_coeff_3p60_vs_3p5`, `q_fig1c_figS1_slope_and_Dmax_anchor`):*
    - `cmp_SI_blue_surface` (`reproduction`, `bundle_geometry_saturation`, $u = 13, \sigma = +1$): **reproduced** (11/11 points inside stroke/calibration bands, though systematically `0.5%–1.5%` below curve centerline).
    - `cmp_SI_green_inside` (`reproduction`, `bundle_geometry_saturation`, $u = 13, \sigma = -1$): **reproduced** (12/12 via input-rounding propagation), **but reveals a hidden equation-vs-plot discrepancy**: combining the digitized blue and green curves across $n = 2..80$ gives $u_{\text{eff}} = 6/(R_+ + R_-) = 12.98 \pm 0.02$ and surface coefficient **$c_{\text{surf}} = 3.605 \pm 0.008$** (independently confirmed by the green curve crossing $R^* = 0.00508$ at $n = 1.512$ and $R^* = 0.0100$ at $n = 1.583$) rather than the printed `3.5` on SI p. 5. Raw `bundle_geometry_saturation` overpredicts the green curve by `+60%` at $n = 1.6$ (`0.0179` vs. `0.0112`), `+16%` at $n = 2$ (`0.0404` vs. `0.0348`), `+6.6%` at $n = 3$, and `+4.4%` at $n = 4$; the automated verdict becomes `reproduced` only because discrete `sigma_inside_only = "-1"` carries a decimal-rounding interval `[-1.5, -0.5]`.
    - `cmp_Rstar_scaling_figS1_line` (`reproduction`, `Rstar_filament_scaling`, fitted prefactor `0.003083`): **not reproduced** (1 reproduced, 3 close, 3 not reproduced) — drawn slope `1.563` matches `1/(2 × 0.327) = 1.529` from Fig. 1C rather than the rounded `5/3 = 1.667`.
    - `cmp_Rstar_scaling_figS1_exp` (`validation`, `Rstar_filament_scaling`, fitted prefactor `0.003661`): **not reproduced** (0/8) due to tight marker intervals and vertical saturation at $n \approx 19.8$.
    - One-to-one pairing with Fig. 1C proves all 8 points use exact anchor **$D_{\max} = D\sqrt{20/n} = 28.22 \pm 0.02\text{ nm}$** (0.07% constancy), settling `q_scaling_approx_provenance`.

- **Fig. S2 (Left: $-13/6$, Right: $-28/13$; p. 7 / SI 3) — Left out of `figures.yaml`:**
  - *Reason:* Radial helical-projection schematic diagrams with no calibrated Cartesian axes; all 12 explicitly printed angular mismatch values (`0°, 4.6°, 9.2°, 13.8°, 9.2°, 4.6°` for $-13/6$ and `6.4°, 10.7°, 2.1°, 6.4°, 2.1°, 10.7°` for $-28/13$) are already recorded directly as `TableData` in `data.yaml` (`mismatch_profile_13_6_figS2`, `mismatch_profile_28_13_figS2`) and reproduced in `cmp_profile_13_6_figS2` and `cmp_profile_28_13_figS2`.
