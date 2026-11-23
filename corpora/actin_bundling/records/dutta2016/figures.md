# Figure Readings & Independent Verification — `dutta2016`

Independent page-by-page inspection and overlay verification of all figures in Dutta, Benetatos & Jho, *EPL* **114**, 28001 (2016) (`corpus/actin_bundling/2016_Dutta.pdf`, PDF pages 2–7, journal pages `28001-p1`–`28001-p6`; PDF page 1 is the IOP publisher cover sheet), verified on 2026-11-23.

All six article pages (PDF pp. 2–7) were rendered at 300 dpi with pixel rulers (`figwork/p2.png`–`figwork/p7.png`) via `boxgraph figure page`. Every quantitative panel in the paper (Figs. 2(a), 2(b), 3(a), 3(b), 3(b) inset, and 4 — 11 curves total, all on PDF p. 5) is already represented by a calibrated 1200–1800 dpi `mode: code` reading in `figures.yaml`; no quantitative panel, inset, curve, or marker series was missed. Every reading was verified via `boxgraph figure check`, full-resolution overlay PNGs (`figwork/<reading_id>.png`), and 3× zoomed cross contact sheets (`figwork/zoom_fig2a.png`, `figwork/zoom_fig2b.png`, `figwork/zoom_fig3_4.png`, `figwork/zoom_updated_4.png`). Four scripts with genuine pixel-level cross misalignments (`fig2b_GA13.py`, `fig3a_bundle_size_vs_GA.py`, `fig3b_bundle_size_vs_lR.py`, and `fig4_finite_boundary.py`) were fixed from pure image geometry and re-verified; the remaining seven accurate scripts were kept unchanged.

- **PDF p. 2 (`28001-p1`, `figwork/p2.png`) — No figures:** Title, authors, abstract, PACS, and Introduction prose only.

- **PDF p. 3 (`28001-p2`, `figwork/p3.png`) — No figures:** Model description, field-theory derivation, and Eqs. (1)–(6) only.

- **PDF p. 4 (`28001-p3`, `figwork/p4.png`, `figwork/fig1_inspect.png`) — Fig. 1(a) & 1(b) — Left out:**
  - *Visual contents:* 3D schematic drawing of charged polymers aligned along the $z$-axis transitioning from (a) disordered in-plane nematic liquid (8 undulating blue chains annotated with symbols $\vec{R}$, $s$, $\vec{r}(s)$, $\mathrm{d}\vec{r}(s)/\mathrm{d}s$, and contour length $L$) to (b) bundled state (3 bundles of 3–4 chains each).
  - *Reason left out:* Purely schematic illustration with no calibrated axes, tick marks, numeric labels, or plotted data.

- **PDF p. 5 (`28001-p4`, `figwork/p5.png`) — Fig. 2(a) — Verified (`fig2a_GA4`, `fig2a_GA13`, `fig2a_GA31`; unchanged):**
  - *Curves:* Dimensionless 2D Fourier potential $\beta V(\mathbf{k})$ vs. dimensionless wave number $k$ at $\Gamma_R = 5$, $\lambda_R = 0.5$, $\lambda_A = 0.3$, $\rho^* = 1$:
    - `fig2a_GA4`: $\Gamma_A = 4$ (red dotted line with circles, no bundles), 11 points $k \in [0.1, 1.4]$.
    - `fig2a_GA13`: $\Gamma_A = 13$ (blue dashed line with triangles, finite bundles), 7 points $k \in [0.15, 1.25]$.
    - `fig2a_GA31`: $\Gamma_A = 31$ (black solid line with squares, macroscopic phase separation), 12 points $k \in [0.1, 1.6]$.
  - *Resolution & calibration:* **1200 dpi**, crop `[1920, 1360, 4920, 3560]`, plot box `2449.5 × 1796.0 px`. Left border ($0.0$) + 3 top x-ticks ($0.5, 1.0, 1.5$): **`x_res = 0.30 px`**; right y-ticks ($4, -4$) + interior horizontal $\beta V = 0$ line: **`y_res = 0.00 px`** (all four borders independently cross-checked).
  - *Overlay verification:* Inspected `figwork/fig2a_GA4.png`, `figwork/fig2a_GA13.png`, `figwork/fig2a_GA31.png`, and `figwork/zoom_fig2a.png`. All 11/11 crosses on `fig2a_GA4` (circle centers), 7/7 on `fig2a_GA13` (dashed centerline, within `<= 1.2 px` of 2D dash-centroid quadratic fit), and 12/12 on `fig2a_GA31` (square centers) sit cleanly on target; all 3 scripts left unchanged.
  - *Reproduction verdicts:* Consistent angular-wavenumber `cmp_fig2a_GA4`, `cmp_fig2a_GA13`, `cmp_fig2a_GA31` **not reproduced**; cyclic/mixed-$2\pi$ variants `cmp_fig2a_GA4_var`, `cmp_fig2a_GA13_var`, `cmp_fig2a_GA31_var` **reproduced**.

- **PDF p. 5 (`28001-p4`, `figwork/p5.png`) — Fig. 2(b) — Verified (`fig2b_GA4`, `fig2b_GA31` unchanged; `fig2b_GA13` fixed):**
  - *Curves:* Stability function $F(\mathbf{k})$ vs. dimensionless wave number $k$ at $\Gamma_R = 5$, $\lambda_R = 0.5$, $\lambda_A = 0.3$, $\rho^* = 1$ (with horizontal threshold line $F = 1$ separating "unstable" and "stable"):
    - `fig2b_GA4`: $\Gamma_A = 4$ (red dotted line with circles), 7 points $k \in [0.1, 1.0]$.
    - `fig2b_GA13`: $\Gamma_A = 13$ (blue dashed line with triangles), 6 points $k \in [0.15, 1.05]$.
    - `fig2b_GA31`: $\Gamma_A = 31$ (black solid line with squares), 11 points $k \in [0.1, 1.4]$.
  - *Resolution & calibration:* **1200 dpi**, crop `[4920, 1360, 8000, 3560]`, plot box `2552.0 × 1796.0 px`. Left/right borders ($0.0, 2.0$) + 2 top x-ticks ($0.5, 1.5$): **`x_res = 0.00 px`**; left y-ticks ($6, 4, 2$), $F = 1$ line, and bottom border ($0$): **`y_res = 0.40 px`**.
  - *Overlay verification & fix:* Inspected `figwork/fig2b_GA4.png`, `figwork/fig2b_GA13.png`, `figwork/fig2b_GA31.png`, `figwork/zoom_fig2b.png`, and `figwork/zoom_updated_4.png`. All 7/7 crosses on `fig2b_GA4` and 11/11 on `fig2b_GA31` sit dead-center on their markers (kept unchanged). In `fig2b_GA13`, the initial column-wise blue-pixel fit suffered a `~4.6 px` downward bias at $k = 0.65$ (`y = 0.9804` vs. `0.9984`) because the black horizontal $F = 1$ line partially occludes the blue dashes near the peak, plus beveled-dash sawtooth (`fit_res = 13.18 px`, `4.9 px` shift at $k = 0.15$). Updated `figures/fig2b_GA13.py` to fit local quadratics to complete unoccluded 2D dash centroids (`area 380–550 px`) plus the bounding-box midpoint of the dash straddling $F = 1$, reducing `fit_res` to `0.01–1.22 px` with all 6/6 crosses centered cleanly on the dashed curve.
  - *Reproduction verdicts:* Consistent `cmp_fig2b_GA4`, `cmp_fig2b_GA13`, `cmp_fig2b_GA31` **not reproduced**; mixed-$2\pi$ variants `cmp_fig2b_GA4_var`, `cmp_fig2b_GA13_var`, `cmp_fig2b_GA31_var` **reproduced**.

- **PDF p. 5 (`28001-p4`, `figwork/p5.png`) — Fig. 3(a) — Verified (`fig3a_bundle_size_vs_GA` fixed):**
  - *Curve:* Dimensionless bundle size vs. attractive interaction strength $\Gamma_A$ at $\Gamma_R = 5$, $\lambda_R = 0.5$, $\lambda_A = 0.3$, $\rho^* = 1$ (solid red piecewise-linear polyline), 12 points $\Gamma_A \in [14, 31]$.
  - *Resolution & calibration:* **1200 dpi**, crop `[2080, 4880, 4880, 6960]`, plot box `2304.0 × 1644.0 px`. Four top x-ticks ($15, 20, 25, 30$): **`x_res = 0.00 px`**; 2 right y-ticks ($20, 10$) + bottom border ($0$): **`y_res = 0.00 px`**.
  - *Overlay verification & fix:* Inspected `figwork/fig3a_bundle_size_vs_GA.png`, `figwork/zoom_fig3_4.png`, and `figwork/zoom_updated_4.png`. Points $\Gamma_A = 14\text{--}29$ sat dead-center, but vertical column slicing caused a `10.0 px` upward miter-join bias at the sharp polyline kink $\Gamma_A = 30$ (`13.54` vs. `13.39`) and a `12.0 px` downward end-cap truncation bias (plus silent `slope_px = 0.0` fallback) at the right endpoint $\Gamma_A = 31$ (`27.60` vs. `27.82`). Updated `figures/fig3a_bundle_size_vs_GA.py` to use horizontal-scan centerline interpolation on steep segments (`span > 16 px`) and clean interior linear centerline extrapolation at the right endpoint; all 12/12 crosses sit dead-center on the polyline and its vertices.
  - *Reproduction verdicts:* Consistent `cmp_fig3a` **not reproduced**; mixed-$2\pi$ variant `cmp_fig3a_var` **reproduced**.

- **PDF p. 5 (`28001-p4`, `figwork/p5.png`) — Fig. 3(b) main panel — Verified (`fig3b_bundle_size_vs_lR` fixed):**
  - *Curve:* Dimensionless bundle size vs. repulsive screening length $\lambda_R$ at $\Gamma_R = 5$, $\Gamma_A = 20$, $\lambda_A = 0.3$, $\rho^* = 1$ (solid red curve, inset region masked out), 10 points $\lambda_R \in [0.40, 0.95]$.
  - *Resolution & calibration:* **1200 dpi**, crop `[4880, 4880, 7760, 6960]`, plot box `2303.0 × 1647.0 px`. Three bottom x-ticks ($0.5, 0.7, 0.9$): **`x_res = 0.00 px`**; 5 left y-ticks ($12, 10, 8, 6, 4$): **`y_res = 0.40 px`**.
  - *Overlay verification & fix:* Inspected `figwork/fig3b_bundle_size_vs_lR.png`, `figwork/zoom_fig3_4.png`, and `figwork/zoom_updated_4.png`. Points $\lambda_R = 0.45\text{--}0.95$ sat dead-center, but on the nearly vertical curved left branch at $\lambda_R = 0.40$, vertical column slicing spanned `135 px` and pulled the cross `38.0 px` too high onto the right edge of the red stroke (`8.251` vs. `7.974`, with `2.1 px` shift at $\lambda_R = 0.42$). Updated `figures/fig3b_bundle_size_vs_lR.py` to switch to horizontal-scan centerline interpolation whenever vertical slice height exceeds `16 px`; all 10/10 crosses now sit dead-center on the stroke.
  - *Reproduction verdicts:* Consistent `cmp_fig3b` **not reproduced**; mixed-$2\pi$ variant `cmp_fig3b_var` **reproduced**.

- **PDF p. 5 (`28001-p4`, `figwork/p5.png`) — Fig. 3(b) inset — Verified (`fig3b_inset_ratio_vs_lR`; unchanged):**
  - *Curve:* Ratio $\text{Bundle size} / \lambda_R$ vs. $\lambda_R$ at $\Gamma_R = 5$, $\Gamma_A = 20$, $\lambda_A = 0.3$, $\rho^* = 1$ (solid red curve), 9 points $\lambda_R \in [0.42, 0.95]$.
  - *Resolution & calibration:* **1800 dpi**, crop `[9000, 7620, 11280, 9300]`, plot box `1634.5 × 1225.5 px`. Three bottom x-ticks ($0.5, 0.7, 0.9$): **`x_res = 0.00 px`**; 3 left y-ticks ($30, 20, 10$): **`y_res = 0.00 px`**.
  - *Overlay verification:* Inspected `figwork/fig3b_inset_ratio_vs_lR.png` and `figwork/zoom_fig3_4.png`. All 9/9 magenta crosses sit cleanly centered on the red curve (`<= 1.4 px` difference from horizontal scan at $\lambda_R = 0.42$); script kept unchanged.
  - *Reproduction verdicts:* Consistent `cmp_fig3b_inset` **not reproduced**; mixed-$2\pi$ variant `cmp_fig3b_inset_var` **reproduced**.

- **PDF p. 5 (`28001-p4`, `figwork/p5.png`) — Fig. 4 — Verified (`fig4_infinite_boundary` unchanged; `fig4_finite_boundary` fixed):**
  - *Curves:* Stability diagram $\lambda_R$ vs. $\Gamma_A$ at $\Gamma_R = 10$, $\lambda_A = 0.1$, $\rho^* = 1$ showing "Nematic liquid", "Finite bundles", and "Infinite bundles" phases:
    - `fig4_finite_boundary`: Nematic liquid / Finite bundles spinodal boundary $F(k) = 1$ at finite $k$ (red dashed line with filled circles), 10 points $\Gamma_A \in [64, 81]$ (the marker at $\Gamma_A = 82$ touches the top border and is excluded).
    - `fig4_infinite_boundary`: Macrophase separation ($F(k=0) = 1$) infinite-bundle boundary (blue solid line with downward triangles), 12 points $\Gamma_A \in [56, 100]$ sampled cleanly on the pure solid blue line between triangle markers.
  - *Resolution & calibration:* **1200 dpi**, crop `[5560, 8080, 8600, 10320]`, plot box `2476.0 × 1793.0 px`. Five bottom x-ticks ($60, 70, 80, 90, 100$): **`x_res = 0.40 px`**; 2 right y-ticks ($0.4, 0.2$) + bottom border ($0.0$): **`y_res = 0.00 px`**.
  - *Overlay verification & fix:* Inspected `figwork/fig4_finite_boundary.png`, `figwork/fig4_infinite_boundary.png`, `figwork/zoom_fig3_4.png`, and `figwork/zoom_updated_4.png`. All 12/12 crosses on `fig4_infinite_boundary` sit dead-center on the solid blue line (kept unchanged). In `fig4_finite_boundary`, per-column medians were pulled upward by steep red dashes directly above the filled circles at $\Gamma_A = 80$ (`5.4 px` high, `0.3303` vs. true circle center `0.3288`) and $\Gamma_A = 81$ (`8.1 px` high, `0.3901` vs. true circle center `0.3879`, visibly sitting near the upper rim of the circle), and $\Gamma_A = 68$ spuriously reported `slope = 0.00`. Updated `figures/fig4_finite_boundary.py` to isolate filled red circle interiors via Euclidean distance transform (`distance_transform_edt(red) >= 10 px`), placing all 10/10 crosses within `<= 0.2 px` of the true circle centers with adjacent-marker local slopes.
  - *Reproduction verdicts:* Analytical $k=0$ boundary `cmp_fig4_infinite` **reproduced**; consistent finite boundary `cmp_fig4_finite` **not reproduced**; mixed-$2\pi$ finite boundary `cmp_fig4_finite_var` **reproduced**.

- **PDF p. 6 (`28001-p5`, `figwork/p6.png`) — No figures:** Salt-screening discussion, Conclusion and discussion, Acknowledgments, and References [1]–[23] only.

- **PDF p. 7 (`28001-p6`, `figwork/p7.png`) — No figures:** References [24]–[45] only.
