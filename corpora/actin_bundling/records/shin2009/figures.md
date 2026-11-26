# Figure Readings — `shin2009`

Calibrated `mode: code` digitization of every quantitative plot across all 4 PDF pages of Shin, Purdy Drew, Bartles, Wong & Grason, *Phys. Rev. Lett.* **103**, 238102 (2009) (`corpus/actin_bundling/2009_Shin.pdf`, pp. 1–4). Fig. 1(c–e) panels were rendered at **1500 dpi** (plot boxes 2086–2963 × 1680–1740 px) and Fig. 3(b) at **1800 dpi** (plot box 3858 × 1570 px). All axes were calibrated directly from detected frame lines and tick marks via `figures/_calib.py`, and all **13 overlay PNGs** in `figwork/` were generated with `boxgraph figure check` and visually verified via contact sheets (`figwork/verify_fig1c.png`, `figwork/verify_fig1d.png`, `figwork/verify_fig1e_fig3b.png`) to confirm that all **206 digitized points** and error-bar caps sit cleanly on their target curves and markers. Saved readings were written to `readings.json`.

- **Page 1 (`238102-1`) — No figures:** Title, abstract, introduction, and experimental methods text only.

- **Fig. 1(a) (p. 2 / `238102-2`) — Left out:**
  - *Reason:* 2D SAXS detector diffraction image of fascin-actin bundles ($R = 0.3$) with schematic $q_z, q_r$ orientation arrows only; no numerical axis ticks or intensity color-bar scale.

- **Fig. 1(b) (p. 2 / `238102-2`) — Left out:**
  - *Reason:* 2D SAXS detector diffraction image of espin-actin bundles ($R = 0.08$) with schematic $q_z, q_r$ orientation arrows only; no numerical axis ticks or intensity color-bar scale.

- **Fig. 1(c) (p. 2 / `238102-2`) — Read (`fig1c_fascin_R0p33`, `fig1c_fascin_R0p2`, `fig1c_fascin_R0p1`, `fig1c_fascin_R0p05`):**
  - *Series (`origin: measured`, 78 points total):* Circularly averaged pseudo-Voigt-background-subtracted SAXS intensity $I(q)$ [a.u.] vs. $q$ [$\text{\AA}^{-1}$] for fascin:actin bundles at 4 molar ratios:
    1. `fig1c_fascin_R0p33` (`I_saxs_fig1c_fascin_R0p33`, magenta curve, $R = 0.33$): 20 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, peaking at $I = 2.852$ ($q = 0.120$) and secondary bump $I = 1.582$ ($q = 0.1335$–$0.1345$).
    2. `fig1c_fascin_R0p2` (`I_saxs_fig1c_fascin_R0p2`, black curve, $R = 0.2$): 20 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, isolated from text labels, arrows, and vertical dashed lines via left-to-right stroke tracking (`track_black_curve`), peaking at $I = 1.880$ ($q = 0.120$).
    3. `fig1c_fascin_R0p1` (`I_saxs_fig1c_fascin_R0p1`, red curve, $R = 0.1$): 19 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, peaking near $q = 0.119$ ($I = 1.322$) with wider vertical intervals reflecting the plotted curve noise.
    4. `fig1c_fascin_R0p05` (`I_saxs_fig1c_fascin_R0p05`, blue curve, $R = 0.05$): 19 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, showing the broad unbundled plateau $I \approx 0.50$–$0.52$ across $q = 0.114$–$0.125\,\text{\AA}^{-1}$.
  - *Calibration (1500 dpi):* Linear x-axis from 6 matched top/bottom inward ticks (`0.11..0.135`) plus right border (`0.140`), max residual `0.43 px` across 2963 px (`0.015%`); linear y-axis from top border (`3.0`), 2 matched left/right ticks (`2.0, 1.0`), and bottom border (`0.0`), max residual `0.17 px` across 1740 px (`0.010%`).
  - *Comparisons & findings (`q_fig1cd_split_axis_and_layer2_peak_offset`):* No runnable SAXS intensity model exists in `shin2009` (`saxs_4sphere_twist_fitting` is non-runnable). Processing model `saxs_layer_line_peak_finder` on `I_saxs_fig1c_fascin_R0p33` yields:
    - `cmp_q1_fig1c_R0p33` (`consistency`): **reproduced** (`0.12041 A^-1` vs. printed `0.120 A^-1`, diff `+0.00041`).
    - `cmp_q2_fig1c_R0p33` (`consistency`): **reproduced** (`0.13393 A^-1` vs. printed `0.134 A^-1`, diff `-0.00007`), even though the drawn vertical dashed guide line sits at `0.1350 A^-1`.

- **Fig. 1(d) (p. 2 / `238102-2`) — Read (`fig1d_espin_R0p2`, `fig1d_espin_R0p05`, `fig1d_espin_R0p02`, `fig1d_pure_Factin`):**
  - *Series (`origin: measured`, `from_source: ref-7`, 78 points total):* Circularly averaged SAXS intensity $I(q)$ [a.u.] vs. $q$ [$\text{\AA}^{-1}$] for espin:actin bundles [7] plotted across a **split linear y-axis** separated by `//` axis breaks:
    1. `fig1d_espin_R0p2` (`I_saxs_fig1d_espin_R0p2`, black curve, $R = 0.2$, upper panel `[100, 200]`): 20 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, sharp bundled peaks at $q = 0.120$ ($I = 185.2$) and $q = 0.135$ ($I = 120.8$).
    2. `fig1d_espin_R0p05` (`I_saxs_fig1d_espin_R0p05`, blue curve, $R = 0.05$, upper panel `[100, 200]`): 20 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, peaks at $q = 0.120$ ($I = 113.6$) and $q = 0.135$ ($I = 74.5$).
    3. `fig1d_espin_R0p02` (`I_saxs_fig1d_espin_R0p02`, brown curve, $R = 0.02$, lower panel `[40, 45]`): 19 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, showing coexistence of broad unbundled shoulder (`0.114–0.116`, $I \approx 46.4$) with emerging sharp bundled peak (`0.120`, $I = 47.48$).
    4. `fig1d_pure_Factin` (`I_saxs_fig1d_pure_Factin`, green curve, pure F-actin, lower panel `[40, 45]`): 19 points on $q \in [0.107, 0.139]\,\text{\AA}^{-1}$, broad native 13/6 maximum at $q = 0.115$ ($I = 43.86$).
  - *Calibration (1500 dpi):* Upper sub-panel x-axis from 6 top ticks (`0.11..0.135`, residual `0.29 px`), lower sub-panel x-axis from 6 bottom ticks (`0.11..0.135`, residual `0.48 px`; top and bottom grids are offset by `12.8 px = 0.00016 A^-1` across the `//` break). Upper y-axis from top border (`200`) and long left tick (`100`); lower y-axis from long left tick (`45`) and bottom border (`40`); the `6.0–7.5 px` vertical offset between left and right ticks is propagated directly into every point's `[lo, hi]` interval.
  - *Comparisons & findings (`q_fig1cd_split_axis_and_layer2_peak_offset`):*
    - `cmp_q1_fig1d_R0p2` (`consistency`): **reproduced** (`0.12005 A^-1` vs. printed `0.120 A^-1`, diff `+0.00005`).
    - `cmp_q2_fig1d_R0p2` (`consistency`): **close** (`0.13477 A^-1` vs. printed `0.134 A^-1`, diff `+0.00077`; vertical dashed guide line also sits at `0.1349 A^-1`).
    - `cmp_q1_fig1d_pure` (`consistency`): **close** (`0.11495 A^-1` vs. printed `0.114 A^-1`, diff `+0.00095`; plotted arrow sits at `0.115 A^-1`).

- **Fig. 1(e) (p. 2 / `238102-2`) — Read (`fig1e_espin`, `fig1e_fascin`):**
  - *Series (`origin: measured`, `error_bars: unstated`, 12 markers total):* Measured F-actin overtwist [degrees] vs. cross-linker molar ratio $R$:
    1. `fig1e_espin` (`overtwist_fig1e_espin`, 7 black filled squares — including cleanly decomposed touching pair at $R = 0.066, 0.077$):
       - $R = 0.0206$: `0.900°` (no protruding error bar; marker extent `[0.868°, 0.933°]`)
       - $R = 0.0500$: `0.948°` (error bar `[0.897°, 1.002°]`)
       - $R = 0.0662$: `0.900°` (marker extent `[0.868°, 0.933°]`)
       - $R = 0.0772$: `0.900°` (marker extent `[0.868°, 0.933°]`)
       - $R = 0.1011$: `0.900°` (error bar `[0.793°, 1.002°]`)
       - $R = 0.1522$: `0.900°` (marker extent `[0.868°, 0.933°]`)
       - $R = 0.3336$: `0.802°` (marker extent `[0.773°, 0.832°]`; upward stem belongs to fascin circle directly above)
    2. `fig1e_fascin` (`overtwist_fig1e_fascin`, 5 red open circles with black vertical error bars):
       - $R = 0.0512$: `0.002°` (no protruding error bar; marker extent `[-0.036°, +0.041°]`)
       - $R = 0.0686$: `0.104°` (error bar `[-0.004°, 0.205°]`, lower cap at `y = 0` axis floor)
       - $R = 0.1012$: `0.302°` (error bar `[0.198°, 0.401°]`)
       - $R = 0.2011$: `0.799°` (error bar `[0.698°, 0.901°]`)
       - $R = 0.3316$: `0.900°` (error bar `[0.798°, 1.002°]`, lower cap occluded by espin square at `0.802°` and mirrored from upper cap)
  - *Calibration (1500 dpi):* Linear x-axis from left border (`0.0`) and 6 unobstructed interior bottom ticks (`0.10, 0.15, 0.20, 0.25, 0.30, 0.35`, excluding the outer-shifted right frame line), residual `0.16 px` across 1818 px (`0.009%`); linear y-axis from top border (`1.2`), 5 left ticks (`1.0..0.2`), and bottom border (`0.0`), residual `3.05 px` across 1680 px (`0.18%`, `0.0022 deg`, due to original raster quantization).
  - *Verdicts & findings (`q_fig1e_espin_point_count_and_high_R_dip`):*
    - `cmp_fig1e_espin_vs_28_13` (`validation`, `overtwist_28_13_plateau` predicting constant `0.9890°`): **not reproduced** overall (`2 reproduced`, `4 close`, `1 not reproduced`): matches both error-barred points ($R = 0.050, 0.101$), lies close (`+0.089°`) to the four un-barred `0.900°` markers, and lies `+0.187°` above the un-barred high-$R$ dip to `0.802°` at $R = 0.334$.
    - `cmp_fig1e_fascin_vs_28_13` (`validation`, `overtwist_28_13_plateau` predicting `0.9890°`): **not reproduced** overall (`1 reproduced`, `1 close`, `3 not reproduced`): quantitatively confirms fascin's continuous low-$R$ transition (`0.002° -> 0.104° -> 0.302°`), approaching the 28/13 plateau near the `6/28 = 0.214` bound fraction (`0.799 ± 0.10°`, `close`) and reproducing it at $R = 0.332$ (`0.900 ± 0.10°`, `reproduced`).

- **Fig. 2(a), 2(b), 2(c) (p. 3 / `238102-3`) — Left out:**
  - *Reason:* Purely schematic illustrations of (a) two actin filaments linked by a cross-link, (b) X,Y spin representation $\hat{\mathbf{S}}_i, \hat{\mathbf{S}}_j$ with filament spacing $D$ and cross-link length $\Delta_{ij}$, and (c) top-view angular layout of the 6 cross-linkers (`0; 4,19; 9,24; 14`) in the 28/13 unit cell; no numeric plot axes.

- **Fig. 3(a) & Inset (p. 3 / `238102-3`) — Left out:**
  - *Reason:* Qualitative mean-field equation-of-state plot ($\rho$ vs. $z$ across $u < u_c$, $u = u_c$, $u > u_c$) and qualitative inset phase diagram ($U$ vs. $z$ with first-order line and critical end point). Neither the main panel nor the inset carries any numeric tick marks or tick labels on either axis, so no calibrated coordinates can be read.

- **Fig. 3(b) (p. 3 / `238102-3`) — Read (`fig3b_xit1000`, `fig3b_xit100`, `fig3b_xit10`):**
  - *Series (`origin: predicted`, 38 points total):* Predicted monomer-monomer overtwist $\operatorname{Im}[\ln g]$ [radians] vs. $(\beta U \rho)^2$:
    1. `fig3b_xit1000` (`overtwist_fig3b_xit1000`, blue top curve, $\xi_t = 1000$): 13 points on $x \in [0.01, 0.25]$ (avoiding exact crossing with $\xi_t = 100$ near $x \approx 0.047$), rising concavely-up from `0.00824 rad` ($x = 0.01$) through `0.07765 rad` ($x = 0.055$) to `0.11372 rad` ($x = 0.25$).
    2. `fig3b_xit100` (`overtwist_fig3b_xit100`, purple middle curve, $\xi_t = 100$): 13 points on $x \in [0.01, 0.25]$, rising from `0.02031 rad` ($x = 0.01$) through `0.07158 rad` ($x = 0.055$) to plateau `0.07834 rad` ($x = 0.25$).
    3. `fig3b_xit10` (`overtwist_fig3b_xit10`, olive bottom curve, $\xi_t = 10$): 12 points on $x \in [0.02, 0.25]$ (below $x = 0.015$ the curve merges into the thick black x-axis arrow), rising from `0.00354 rad` ($x = 0.02$) to plateau `0.00789 rad` ($x = 0.25$).
  - *Calibration (1800 dpi):* Linear x-axis strictly from the 5 major ticks (`0.05..0.25`, residual `0.40 px` across 3085 px = `0.013%`); linear y-axis strictly from the 6 major ticks (`0.02..0.12`, residual `0.30 px` across 1308 px = `0.023%`). Extrapolating the tick grid to `y = 0` places the true plot origin at row `1886.3` (matching quadratic extrapolation of the olive curve `1886.1` to `0.2 px`), whereas the overlaid black x-axis arrow sits `10.8 px = 0.0008 rad` higher at row `1875.5`.
  - *Verdicts & settled/sharpened open question (`q_fig3b_overtwist_unstated_formula_and_discrepancy`, `q_fig3b_digitized_quadratic_abscissa_and_asymptotes`):*
    - Exact connected $\mathcal{O}(U^2)$ perturbation theory (`overtwist_fig3b_perturbative`): **not reproduced** on all 3 curves — `cmp_fig3b_pert_xit1000` (`0/13`: `0.0205` vs. `0.1137 rad` at $x = 0.25$, 5.5× low), `cmp_fig3b_pert_xit100` (`0/13`: `0.8709` vs. `0.0783 rad`, 11× high), `cmp_fig3b_pert_xit10` (`0/12`: `0.1577` vs. `0.00789 rad`, 20× high).
    - Exact converged (`n_modes = 16`, 8-digit converged) 1D Fourier transfer matrix (`overtwist_fig3b_transfer_matrix`): **not reproduced** on all 3 curves — `cmp_fig3b_tm_xit1000` (`0/13`: `0.0610` vs. `0.1137 rad`, opposite ordering to $\xi_t = 100$), `cmp_fig3b_tm_xit100` (`6 close` on $x = 0.11$–$0.25$ where transfer matrix saturates at $\pi/39 = 0.0805\text{ rad}$ vs. plotted `0.0783 rad`, `7 not reproduced` on $x = 0.01$–$0.09$ where it is concave down and 2.1× above plotted `0.0203 rad`), `cmp_fig3b_tm_xit10` (`0/12`: `0.0626` vs. `0.00789 rad`, 8× high).
    - Empirically, all three digitized Fig. 3(b) curves fit $\arctan\!\left(\frac{A x^2}{1 + B x^2}\right)$ in the plotted abscissa $x$ to `0.0003–0.0009 rad` (12–30× better than any $\mathcal{O}(x)$ ratio), with asymptotes `0.00786`, `0.07864` (exact 10.0× ratio), and `0.1163 rad`, confirming that Fig. 3(b) cannot be reproduced by the paper's stated $\mathcal{O}(U^2)$ expansion of $g(\ell_0)$ versus $(\beta U \rho)^2$.

- **Page 4 (`238102-4`) — No figures:** Theory discussion, conclusions, acknowledgments, and references only.
