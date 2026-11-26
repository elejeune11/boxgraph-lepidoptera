# Figure Inventory for `bathe2008`

Complete page-by-page, panel-by-panel figure audit and calibrated digitization of `corpus/actin_bundling/2008_Bathe.pdf` (Bathe, Heussinger, Claessens, Bausch & Frey, *Biophysical Journal* **94**(8), 2955–2964, 2008; PDF pp. 1–10 main article, PDF p. 11 Supplementary Material converted from `mmc1.doc`).

All **55 calibrated figure readings** (**348 digitized points** across all 8 quantitative plot panels) use fast `mode: code` Python readers under [`figures/`](figures/) delegating to [`figures/_helpers.py`](figures/_helpers.py) (which operate directly on the rendered image array with no PDF re-rendering). Every single reading was checked via `.venv/bin/boxgraph figure check` into `figwork/overlays/<id>.png`, visually verified on compact contact sheets (`figwork/sheet_1.png`–`sheet_7.png`, `figwork/verify_new_1.png`–`verify_new_3.png`, and `figwork/verify_fixed_zoom.png`) with zero defects, wired into `data.yaml`, `runs.yaml`, and `comparisons.yaml`, and saved to [`readings.json`](readings.json) (55/55 fresh matching recipes).

## Summary of Dedicated `digitize:bathe2008` Pass Additions & Fixes

- **Added 14 missing quantitative readings (+94 points, bringing 41 → 55 readings and 254 → 348 points)**:
  - **Fig. 3A**: Added `fig3a_a1e_1` ($\alpha = 10^{-1}$, 6 points at $N \in \{49, 100, 150, 225, 300, 360\}$ where the thin curve cleanly separates above the thick decoupled $\kappa_B^* = N$ line).
  - **Fig. 3B**: Added all 6 remaining solid FE curves `fig3b_n9`, `fig3b_n16`, `fig3b_n36`, `fig3b_n49`, `fig3b_n64`, `fig3b_n81` (10 points each = 60 points), completing all 9/9 solid FE curves ($N \in \{4, 9, 16, 25, 36, 49, 64, 81, 100\}$) alongside all 3/3 dotted Timoshenko curves.
  - **Fig. 3B inset**: Added the two degenerate $N = 4$ crossover markers `fig3b_inset_n4_pinned` (open square, $\alpha^{**} = 2.981\;[2.656, 3.346]$ vs Eq. 2 value $3.0$) and `fig3b_inset_n4_clamped` (open circle, $\alpha^{**} = 11.591\;[10.326, 13.010]$ vs Eq. 2 value $12.0$), directly confirming `q_fig3b_inset_n4_degenerate`.
  - **Fig. 4**: Preserved pooled `fig4_exp_points` (all 20 measured symbols) and added the 5 concentration-/spacing-resolved measured symbol series (`fig4_exp_d40` 5 squares, `fig4_exp_d56` 3 circles, `fig4_exp_d68` 4 diamonds, `fig4_exp_d225` 4 up-triangles, `fig4_exp_d412` 4 down-triangles) classified purely from marker geometry, recording each series' `crosslink_spacing` condition (`40 nm`, `56 nm`, `68 nm`, `225 nm`, `412 nm`) so future cross-paper validations can recover both $\delta$ and $L = \sqrt{(L^2/\delta)\delta}$.
- **Fixed calibration/reading defects & performance**:
  - **Fig. S2B (`figS2b_m4`) hollow-marker split bug**: Cross-checking identical $c$ values plotted in Fig. S2A vs Fig. S2B revealed that `_clusters(..., gap=12)` split the ~53 px hollow red octagon markers at $x = 16$ ($\phi = 2$) and $x = 6400$ ($\phi = 40$) into top/bottom rims, shifting readings by $+6.35\%$ and $-8.10\%$. Narrowing the centered column slice to `ci ± 12 px` with `gap = 70 px` centers every hollow marker cleanly, bringing `figS2b_m4` agreement with `figS2a_m4` to $0.72\%$–$1.32\%$ (and all 25 S2A/S2B marker pairs to $\le 1.32\%$).
  - **Fig. 3A tick-mark margin**: Increased vertical border exclusion in `read_fig3a_curve` from `25 px` to `50 px` to clear the `32 px` inward major x-tick marks at $N = 4, 9, 100, 300$.
  - **Reader speedup (~20–40×)**: Replaced full-image `float64` `image.mean(axis=2)` conversions, 2D `binary_opening`, and full-frame `distance_transform_edt` in `figures/_helpers.py` with `uint16`/`uint8` masks, separable 1D box morphology, and per-component padded-bbox EDT (`read_fig4_exp_markers` dropped from `28.6 s` to `0.78 s`).

---

## Panel-by-Panel Inventory (Pages 1–11)

### Figures digitized (in `figures.yaml`, `figures/`, and `readings.json`)

- **Fig. 3A** (PDF p. 4, printed p. 2958): Theoretical normalized bending stiffness $\kappa_B^*$ vs filament number $N$ at fixed $\alpha$.
  - *Readings (5 readings, 38 points)*: `fig3a_a1e_1` ($\alpha = 10^{-1}$, 6 pts at $N \in \{49, 100, 150, 225, 300, 360\}$), `fig3a_a1e0` ($\alpha = 10^0$, 8 pts), `fig3a_a1e1` ($\alpha = 10^1$, 8 pts), `fig3a_a1e2` ($\alpha = 10^2$, 8 pts), and `fig3a_a1e3` ($\alpha = 10^3$, 8 pts).
  - *Calibration*: 1200 dpi vector crop, plot box `2932 × 2196 px`, log10 linear fit to paired inward top/bottom x-ticks ($N \in [1, 400]$) and left/right y-ticks ($\kappa_B^* \in [1, 5 \times 10^5]$); max tick misfit `(2.81, 3.18) px`.
  - *Left out in this panel*: For $\alpha = 10^{-1}$, $N < 49$ is merged into the thick decoupled ($\kappa_B^* = N$) line and omitted; $\alpha = 10^4$ differs from the thick fully coupled ($\kappa_B^* = N^2$) line by $\le 3.6\%$ ($\le 6\text{ px}$, less than half the thick-line stroke width) even at $N = 360$ and remains merged across the entire plot, so it is omitted to avoid ambiguous stroke separation. The two thick analytical boundary lines ($\kappa_B^* = N$ and $\kappa_B^* = N^2$) are exact identities rather than FE data.

- **Fig. 3B** (PDF p. 4, printed p. 2958): Theoretical normalized bending stiffness $\kappa_B^*$ vs fiber coupling parameter $\alpha$ at constant $N \in \{4, 9, 16, 25, 36, 49, 64, 81, 100\}$ plus dotted Timoshenko curves.
  - *Readings (12 readings, 100 points)*: All 9 solid FE curves `fig3b_n4`, `fig3b_n9`, `fig3b_n16`, `fig3b_n25`, `fig3b_n36`, `fig3b_n49`, `fig3b_n64`, `fig3b_n81`, `fig3b_n100` (10 pts each; for `fig3b_n100`, $\alpha = 3$ is skipped because it touches the bottom-right corner of the inset border and $\alpha = 5$ is sampled instead) plus all 3 dotted Timoshenko curves `fig3b_timo_n4` (2 pts), `fig3b_timo_n25` (3 pts), `fig3b_timo_n100` (5 pts).
  - *Calibration*: 1200 dpi vector crop, plot box `2932 × 2200 px`, log10 linear fit to paired inward x-ticks ($\alpha \in [10^{-2}, 4 \times 10^3]$) and y-ticks ($\kappa_B^* \in [1, 2 \times 10^4]$); max tick misfit `(2.27, 2.99) px`.

- **Fig. 3B inset** (PDF p. 4, printed p. 2958): Factor-of-2 regime crossover values $\alpha^{**}$ vs $N$ for pinned (open squares) and clamped (open circles) boundary conditions.
  - *Readings (6 readings, 34 points)*: Separated branch sweeps for $N \in \{9, 16, 25, 36, 49, 64, 81, 100\}$ (`fig3b_inset_pinned_dec`, `fig3b_inset_pinned_fc`, `fig3b_inset_clamped_dec`, `fig3b_inset_clamped_fc`; 8 pts each) plus the two degenerate $N = 4$ markers where $2N = N^2/2 = 8$ (`fig3b_inset_n4_pinned`, `fig3b_inset_n4_clamped`; 1 pt each).
  - *Calibration*: 2400 dpi vector crop, plot box `1984 × 1468 px`, log10 linear fit to bottom x-ticks ($N \in [4, 100]$; max misfit `3.39 px`) and the two paired major y-ticks ($10^0, 10^2$; nominal `1.0 px` uncertainty). The four solid straight asymptotic guide lines ($N$-independent and linear-in-$N$ scaling guides) are reference guides rather than independent simulation data.

- **Fig. 4** (PDF p. 6, printed p. 2960): Experimental bending stiffness $\kappa_B$ [$\text{N}\cdot\text{m}^2$] vs $L^2/\delta$ [$\text{m}$] for fascin-cross-linked F-actin bundles ($N = 27 \pm 6$, Claessens et al. 2006, Ref. 29) and Eq. (2) ($c = 5$) theoretical curves.
  - *Readings (9 readings, 64 points)*: All 3 theoretical curves — `fig4_theory_n27` (solid $N = 27$, 8 pts), `fig4_theory_n21` (lower dashed $N = 21$, 8 pts), `fig4_theory_n33` (upper dashed $N = 33$, 8 pts); pooled measured dataset `fig4_exp_points` (all 20 filled black symbols, `origin: measured`, `from_source: ref-29`); and all 5 symbol-resolved measured concentration/spacing series (`fig4_exp_d40`: 5 squares, $\delta = 40\text{ nm}$; `fig4_exp_d56`: 3 circles, $\delta = 56\text{ nm}$; `fig4_exp_d68`: 4 diamonds, $\delta = 68\text{ nm}$; `fig4_exp_d225`: 4 up-triangles, $\delta = 225\text{ nm}$; `fig4_exp_d412`: 4 down-triangles, $\delta = 412\text{ nm}$; no error bars are plotted on symbols in Fig. 4).
  - *Calibration*: 1200 dpi vector crop, plot box `2944 × 2208 px`, log10 linear fit to paired inward x-ticks ($L^2/\delta \in [2 \times 10^{-4}, 10^0]\text{ m}$) and y-ticks ($\kappa_B \in [5 \times 10^{-25}, 5 \times 10^{-23}]\text{ N}\cdot\text{m}^2$); max tick misfit `(2.72, 1.77) px`.

- **Fig. 5** (PDF p. 6, printed p. 2960): Bending stiffness state diagram ($\alpha$ vs $N$) for seven in vivo cytoskeletal bundle systems.
  - *Readings (7 readings, 14 points)*: Lower-left $(N_{\min}, \alpha_{\min})$ and upper-right $(N_{\max}, \alpha_{\max})$ box wall centers for all seven colored boxes: `fig5_box_a` (blue, horseshoe crab sperm acrosome), `fig5_box_b` (red, stereocilia), `fig5_box_c` (orange, brush-border microvilli), `fig5_box_d` (cyan, stress fibers), `fig5_box_e` (green, filopodia), `fig5_box_f` (magenta, *Drosophila* neurosensory bristles), and `fig5_box_g` (yellow, outer pillar cell MT bundles).
  - *Calibration*: 1200 dpi crop of embedded figure, plot box `2631 × 2628 px`, log10 linear fit to paired inward x-ticks ($N \in [10^0, 10^4]$) and y-ticks ($\alpha \in [10^{-4}, 10^4]$); max tick misfit `(3.35, 8.61) px` (propagated into corner intervals via `extra = 2 * y_res`). The two dashed regime boundary lines ($\alpha = 1$ and $\alpha = N$) are exact analytical scaling lines.

- **Fig. S1** (PDF p. 11, Supplementary Material): Unnormalized FE bundle bending stiffness $\kappa_B$ vs $\alpha$ at $M = 8$ ($2M = 16$ layers, $N = 256$) for all six fracture densities $\phi \in \{1, 2, 4, 8, 20, 40\}$.
  - *Readings (6 readings, 48 points)*: `figS1_phi1` (black), `figS1_phi2` (blue), `figS1_phi4` (red), `figS1_phi8` (green), `figS1_phi20` (orange/yellow), `figS1_phi40` (magenta), 8 pts each at $\alpha \in \{10^{-1}, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6\}$.
  - *Calibration*: 1200 dpi vector crop, plot box `2866 × 2150 px` spanning $\alpha \in [10^{-3}, 10^7]$ and $\kappa_B \in [10^1, 10^4]$; independently verified against all 89 bottom x-ticks and 28 right y-ticks with max tick residual `(0.50, 1.31) px`.

- **Fig. S2A** (PDF p. 11, Supplementary Material): Fitted Eq. (2) parameter $c$ vs fracture density $\phi \in \{2, 4, 8, 20, 40\}$ before renormalization for all five bundle sizes $M \in \{2, 4, 8, 16, 32\}$.
  - *Readings (5 readings, 25 points)*: `figS2a_m2` (black squares), `figS2a_m4` (red octagons/circles), `figS2a_m8` (blue diamonds), `figS2a_m16` (green triangles), `figS2a_m32` (cyan crosses), 5 pts each.
  - *Calibration*: 1200 dpi vector crop, plot box `2845 × 2140 px` spanning $\phi \in [1, 50]$ and $c \in [10^{-1}, 5 \times 10^3]$; verified against all 12 x-ticks and 39 y-ticks with max tick residual `(0.96, 1.19) px`.

- **Fig. S2B** (PDF p. 11, Supplementary Material): Fitted Eq. (2) parameter $c$ vs renormalized fracture coordinate $(\phi/\phi^*)^2$ with $\phi^* = M/8$ for all five bundle sizes $M \in \{2, 4, 8, 16, 32\}$.
  - *Readings (5 readings, 25 points)*: `figS2b_m2`, `figS2b_m4`, `figS2b_m8`, `figS2b_m16`, `figS2b_m32`, 5 pts each at $(\phi / (M/8))^2$.
  - *Calibration*: 1200 dpi vector crop, plot box `2902 × 2176 px` spanning $(\phi/\phi^*)^2 \in [10^{-2}, 10^5]$ and $c \in [10^{-1}, 10^4]$; verified against all 62 x-ticks and 45 y-ticks with max tick residual `(0.64, 0.80) px`. All 25 marker $c$ values agree with their Fig. S2A counterparts to $0.00\%$–$1.32\%$.

### Figures not digitized (and why)

- **Fig. 1A–D** (PDF p. 2, printed p. 2956): Electron micrographs of F-actin bundles — (A) bullfrog saccule stereocilia (~60 stereocilia), (B) mouse melanoma filopodium, (C) epithelial microvilli, and (D) *Drosophila* neurosensory bristles. No numeric axes, scale-bar plot data, or quantitative measurements.
- **Fig. 2A–B** (PDF p. 3, printed p. 2957): Theoretical bundle model schematics — (A) cross-section ($N = 16$) and side-view spring schematic ($d_f, \delta, L, k_\times, t$); (B) qualitative deformed backbone and axial displacement $u^{(k)}(x, y)$ profiles in decoupled vs fully coupled bending. Qualitative illustration with no numeric axes.
- **Pages 1, 5, 7, 8, 9, 10**: Text, equations, Appendix derivation, and References only — no figures.
