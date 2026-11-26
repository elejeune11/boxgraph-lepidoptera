# Digitized Figures — Tang & Janmey (1996)

All 8 figures in Tang & Janmey (1996, *J. Biol. Chem.* 271(15), 8556–8563, PDF pages 2–6) comprise 13 quantitative KaleidaGraph plots (50 distinct curve/fit-line readings, 434 total digitized points). The paper contains **no** micrographs, gel images, schematics, or non-numeric panels (pages 1, 7, and 8 contain only text and equations).

## Calibration & Extraction Method

- Every panel crop was rendered at **1200 dpi** (plot-box widths 2746–2951 px, heights 1938–2209 px — well above 1500 px on the shorter axis, with full tick labels intact).
- Each panel was calibrated in image space via an automatic 2D affine transform (`figures/_common.py`, `figures/_panels.py`) fitted to all detected inward major/minor ticks on all four plot borders, correcting the ~0.05°–0.25° monochrome scan skew inherent to the 300-dpi source scans. Maximum tick-fitting residuals across all 13 panels are **1.47–2.98 px (x)** and **1.70–3.63 px (y)** at 1200 dpi (0.37–0.91 native 300-dpi scan pixels).
- Filled/hollow symbols, overlapping marker clusters, and plotted regression lines were extracted purely from the passed image array (`figures/_fig1.py`–`_fig8.py`) using 2×-downsampled morphological/cavity detection, full-resolution geometric symbol templates (`circle`, `square`, `diamond`, `uptri`, `downtri`, `dotcircle`, and filled variants), ink-connectivity graph tracing, and vertical dark-run tracking.
- Every single one of the 50 readings was individually checked with `boxgraph figure check` and visually inspected on compact overlay contact sheets to verify dead-center marker/line alignment.

---

## Per-Panel Summary

- **Fig. 1 (p. 2 / printed p. 8557)** — *90° light scattering vs cation concentration (0.1–100 mM, log–log), 0.5 mg/ml F-actin, 150 mM KCl, pH 7.2.*
  - Calibration: box `2797 × 1985 px`, `resid_x = 2.46 px`, `resid_y = 2.80 px`.
  - Readings (8 curves, 63 points): `fig1_arg6` (●, 9 pts), `fig1_his6` (▲, 8 pts), `fig1_his5` (▽, 7 pts), `fig1_cohex` (◇, 8 pts), `fig1_co2` (▼, 7 pts), `fig1_mn2` (△, 8 pts), `fig1_ca2` (□, 7 pts), `fig1_mg2` (○, 9 pts).
  - Notes & cross-checks: Sharp transitions bracket every threshold stated in the p. 2 text: `Arg6` `0.30 → 0.35 mM` (text `0.3 mM`), `His6` `0.80 → 0.89 mM` (`0.8 mM`), `His5` `1.25 → 1.49 mM` (`1.5 mM`), `[Co(NH3)6]3+` `2.99 → 3.19 mM` (`3 mM`), `Co2+` `4.97 → 5.77 mM` (`5 mM`), `Mn2+` `6.06 → 7.96 mM` (`7 mM`), `Ca2+` `14.93 → 19.87 mM` (`20 mM`), `Mg2+` `28.40 → 31.58 mM` (`30 mM`). `His5` includes a baseline point at `0.6 mM` overlapping `His6`; `Ca2+` starts at `5 mM` right beneath the `Co2+` marker; `Mg2+` jumps directly from `10` to `20 mM` (no `15 mM` circle).

- **Fig. 2a (p. 3 / printed p. 8558)** — *Light scattering vs `[Cu2+]` (0–1 mM, linear–log) for 0.5 mg/ml F-actin with and without 20 mM MgCl2.*
  - Calibration: box `2942 × 2081 px`, `resid_x = 2.43 px`, `resid_y = 1.82 px`.
  - Readings (2 curves, 18 points): `fig2a_cu_with_20mM_mg` (●, 8 pts, onset between `0.025` and `0.05 mM` Cu2+, plateau `1.62`), `fig2a_cu_no_mg` (○, 10 pts, onset between `0.30` and `0.40 mM` Cu2+, plateau `1.11`).

- **Fig. 2b (p. 3 / printed p. 8558)** — *Light scattering vs `[Mg2+]` (0–50 mM, linear–log) for 0.5 mg/ml F-actin preincubated with 0, 10, 50, or 100 uM CuCl2.*
  - Calibration: box `2942 × 2081 px`, `resid_x = 2.28 px`, `resid_y = 2.91 px`.
  - Readings (4 curves, 37 points): `fig2b_mg_no_cu` (▽, 10 pts), `fig2b_mg_10uM_cu` (□, 10 pts), `fig2b_mg_50uM_cu` (○, 9 pts), `fig2b_mg_100uM_cu` (◇, 8 pts).
  - Notes: Preincubation with `50` and `100 uM` CuCl2 lowers the Mg2+ bundling transition from `25–30 mM` (`0` and `10 uM` Cu) to `10–15 mM` (`50 uM`) and `5–10 mM` (`100 uM`), while lowering plateau intensity from `~2.0` (`no Cu`) and `~1.3–1.6` (`10 uM`) to `~0.75–0.83`.

- **Fig. 3a (p. 3 / printed p. 8558)** — *Light scattering vs `[Lys18]` (0–100 uM, linear–log) at 6.0, 2.0, and 0.1 mg/ml F-actin (150 mM KCl).*
  - Calibration: box `2951 × 2209 px`, `resid_x = 2.11 px`, `resid_y = 2.13 px`.
  - Readings (3 curves, 33 points): `fig3a_actin_6p0` (□, 11 pts), `fig3a_actin_2p0` (◇, 12 pts), `fig3a_actin_0p1` (○, 10 pts).
  - Notes: The `0.1 mg/ml` trace consists of three closely overlapping open-circle replicate curves; `fig3a_actin_0p1` records their midpoint and vertical envelope `[lo, hi]` across the rising transition (`16.7–27.8 uM`). The three dotted vertical onset lines lie at `[Lys18] = 20.71 uM` (`0.1 mg/ml`), `33.50 uM` (`2.0 mg/ml`), and `55.69 uM` (`6.0 mg/ml`).

- **Fig. 3b (p. 3 / printed p. 8558)** — *Bundling onset `[Lys18]` (uM) vs `[Actin]` (0–150 uM, linear–linear) with printed linear fit `[Lys18] = 19.4 + 0.26 x [Actin]`.*
  - Calibration: box `2951 × 2201 px`, `resid_x = 2.38 px`, `resid_y = 3.54 px`.
  - Readings (2 series, 12 points): `fig3b_lys18_points` (○, 6 measured points) and `fig3b_lys18_fit_line` (plotted fit line sampled at 6 marker-free actin concentrations `10, 30, 60, 80, 110, 130 uM`).
  - Quantitative findings & model comparisons:
    - Measured points: `(2.31, 19.32)`, `(11.90, 22.21)`, `(23.94, 24.94)`, `(47.71, 33.28)`, `(95.17, 44.43)`, `(142.83, 55.56)` uM. The x-values match `0.1, 0.5, 1.0, 2.0, 4.0, 6.0 mg/ml` divided by `42.0 kDa` (`2.38–142.86 uM`) to ~0.1%, directly confirming `M_actin_added = 42000 g/mol`.
    - OLS on the 6 digitized points gives `[Lys18] = 19.32 + 0.2589 * [Actin]` (`R^2 = 0.9959`, RMS `0.84 uM` vs printed equation).
    - `cmp_fig3b_fit_line_repro`: **`reproduced`** (6/6 points, max `|diff| = 0.18 uM` at `130 uM`, reflecting the `0.2589 → 0.26` slope rounding).
    - `cmp_fig3b_fit_quality`: `max |diff| = 1.48 uM` (`+4.7%` at `47.7 uM` actin), reflecting biological scatter around the regression line.
    - `cmp_eq4_vs_fig3b_points`: Appendix Eq. 4 (`K = 10^3, V_P = 1.2 M^-1, q_net = 11, c_18 = 19.4 uM`, slope `0.2561`) reproduces 4/6 digitized points including `6.0 mg/ml` (`55.97` vs `55.56 uM`, `0.7%`), with largest residual `+1.67 uM` (`+5.3%`) at `47.7 uM`.

- **Fig. 4 (p. 4 / printed p. 8559)** — *Light scattering vs `[Lys18]` (0–120 uM, linear–linear) for 0.2 mg/ml F-actin at six gelsolin:actin molar ratios (`no gel`, `1:500`, `1:100`, `1:50`, `1:30`, `1:10`).*
  - Calibration: box `2797 × 1989 px`, `resid_x = 2.98 px`, `resid_y = 3.49 px`.
  - Readings (6 curves, 55 points): `fig4_no_gel` (△, 9 pts), `fig4_gel_1_500` (○, 9 pts), `fig4_gel_1_100` (◇, 9 pts), `fig4_gel_1_50` (▽, 9 pts), `fig4_gel_1_30` (□, 9 pts), `fig4_gel_1_10` (⊙, 10 pts, includes extra point at `100.1 uM`).
  - Notes: Onset remains near `~22 uM` for `no gel` through `1:50` (`>= 135–140 nm` average filament length) while maximal scattering intensity progressively decreases (`1.00 → 0.93 → 0.78 → 0.78 → 0.72 → 0.50`) and the transition shifts rightward and broadens sharply at `1:30` and `1:10`.

- **Fig. 5a (p. 5 / printed p. 8560)** — *Light scattering vs `[Co(NH3)6^3+]` (0–5 mM, linear–log) for 0.2 mg/ml F-actin at 30, 50, 100, and 150 mM KCl.*
  - Calibration: box `2746 × 1938 px`, `resid_x = 2.90 px`, `resid_y = 1.87 px`.
  - Readings (4 curves, 47 points): `fig5a_kcl_30mM` (○, 10 pts), `fig5a_kcl_50mM` (□, 12 pts), `fig5a_kcl_100mM` (△, 11 pts), `fig5a_kcl_150mM` (◇, 14 pts).
  - Notes: Sharp onset transitions shift monotonically with KCl (`~0.8 mM` at 30 mM, `~1.1–1.2 mM` at 50 mM, `~1.9–2.0 mM` at 100 mM, `~2.9–3.0 mM` at 150 mM KCl), consistent with preparation 1 in Fig. 5b.

- **Fig. 5b (p. 5 / printed p. 8560)** — *Bundling onset `[Mn2+]` and `[Co(NH3)6^3+]` vs `[KCl]` (30–150 mM, log–log) with printed power-law fits `[Mn2+] = 3.11 x [KCl]^0.20` and `[Co(NH3)6^3+] = 0.068 x [KCl]^0.74`.*
  - Calibration: box `2747 × 1941 px`, `resid_x = 2.65 px`, `resid_y = 1.70 px`.
  - Readings (6 series, 26 points): `fig5b_mn_points` (●, 4 pts), `fig5b_mn_fit_line` (5 sampled points at `35, 45, 65, 85, 125 mM`), `fig5b_cohex_prep1_circle` (○, 4 pts), `fig5b_cohex_prep2_uptri` (△, 4 pts), `fig5b_cohex_prep3_downtri` (▽, 4 pts), `fig5b_cohex_fit_line` (5 sampled points).
  - Quantitative findings & model comparisons:
    - `Mn2+` points: `(29.99, 6.01)`, `(50.03, 7.00)`, `(99.85, 7.49)`, `(150.0, 8.50)` mM; OLS gives `3.12 * [KCl]^0.198`. `cmp_fig5b_mn_fit_line_repro`: **`reproduced`** (5/5, max `|diff| = 0.069 mM`); `cmp_fig5b_mn_fit_quality`: **`close`** (max `|diff| = 0.32 mM` at 100 mM KCl). At 150 mM KCl, measured `8.50 mM` confirms the ~21% gap above Fig. 1's `7 mM`.
    - `Co(NH3)6^3+` points across 3 independent preps: `0.80–0.90 mM` (30 mM), `1.10–1.29 mM` (50 mM), `1.99–2.00 mM` (100 mM), `2.70–3.00 mM` (150 mM); pooled 12-pt OLS gives `0.0671 * [KCl]^0.743`. `cmp_fig5b_cohex_fit_line_repro`: **`reproduced`** (5/5, max `|diff| = 0.016 mM`); `cmp_fig5b_cohex_fit_prep1..3` quantify 11–17% prep-to-prep scatter around the pooled power law.
    - Trace Manning Eq. 3 shape tests (`cmp_eq3_mn_fig5b` with `N = 2` and `cmp_eq3_cohex_prep1..3_fig5b` with `N = 3`): **`not reproduced`** across every point (`max |diff| = 25.6 mM` for Mn2+, `12.9–13.7 mM` for CoHex), confirming quantitatively on the full curves that Eq. 3's trace-limit power law is far steeper than measured.

- **Fig. 6a (p. 5 / printed p. 8560)** — *Disassembly of 0.5 mg/ml F-actin + 100 uM Lys18 bundles vs `[Nucleoside Triphosphate]` (0–3 mM, linear–linear) for ATP, CTP, and GTP.*
  - Calibration: box `2944 × 2078 px`, `resid_x = 2.47 px`, `resid_y = 3.63 px`.
  - Readings (3 curves, 30 points): `fig6a_atp` (□, 10 pts), `fig6a_ctp` (▽, 10 pts), `fig6a_gtp` (◇, 10 pts).
  - Notes: All three NTPs disassemble bundles sharply over `1.0–1.8 mM` with nearly identical midpoints (`~1.3 mM` for GTP, `~1.4–1.5 mM` for ATP and CTP), demonstrating charge-driven rather than nucleotide-specific disassembly.

- **Fig. 6b (p. 5 / printed p. 8560)** — *Disassembly of 0.3 mg/ml F-actin + 50 uM Lys18 bundles vs `[Adenosine Phosphate]` (0.1–100 mM, log–log) for ATP, ADP, and AMP.*
  - Calibration: box `2942 × 2078 px`, `resid_x = 1.96 px`, `resid_y = 2.40 px`.
  - Readings (3 curves, 27 points): `fig6b_atp` (○, 8 pts), `fig6b_adp` (◇, 6 pts), `fig6b_amp` (△, 13 pts).
  - Notes: Disassembly midpoints scale steeply with anion valence: `~2.0 mM` for ATP, `~5.0 mM` for ADP, and `~20–25 mM` for AMP.

- **Fig. 7 (p. 5 / printed p. 8560)** — *Light scattering vs polycationic peptide concentration (0.0001–0.1 mM, log–log) for 0.1 mg/ml F-actin (150 mM KCl): Lys42, MARCKS peptide, and Lys18.*
  - Calibration: box `2801 × 1989 px`, `resid_x = 2.03 px`, `resid_y = 1.79 px`.
  - Readings (3 curves, 28 points): `fig7_lys42` (○, 11 pts, `0.0001–0.010 mM`), `fig7_marcks` (◇, 9 pts, `0.00067–0.012 mM`), `fig7_lys18` (□, 8 pts, `0.0055–0.055 mM`).
  - Notes: Bundling onset shifts down roughly an order of magnitude per valence step: `~10–20 uM` (`Lys18`), `~1.3–2 uM` (`MARCKS`), and `< 0.2 uM` (`Lys42`, rising steadily from `0.1 uM` to saturation at `1.5 uM`).

- **Fig. 8a (p. 6 / printed p. 8561)** — *Light scattering vs `[Calponin]` (0–6 uM, linear–log) for 2.5 uM F-actin at 30, 50, and 75 mM KCl.*
  - Calibration: box `2942 × 2080 px`, `resid_x = 2.30 px`, `resid_y = 2.43 px`.
  - Readings (3 curves, 29 points): `fig8a_kcl_30mM` (○, 8 pts), `fig8a_kcl_50mM` (◇, 8 pts), `fig8a_kcl_75mM` (▽, 13 pts).
  - Notes: Increasing KCl shifts the calponin bundling transition from `~1.2–1.8 uM` (30 mM KCl) to `~1.8–2.3 uM` (50 mM KCl) and `~3.2–4.1 uM` (75 mM KCl).

- **Fig. 8b (p. 6 / printed p. 8561)** — *Disassembly of 2.5 uM F-actin + 3 uM calponin bundles (50 mM KCl) vs `[Adenosine Phosphate]` (0.5–21 mM, log–log) for AMP, ADP, and ATP.*
  - Calibration: box `2942 × 2080 px`, `resid_x = 1.47 px`, `resid_y = 1.92 px`.
  - Readings (3 curves, 29 points): `fig8b_amp` (○, 10 pts), `fig8b_adp` (◇, 10 pts), `fig8b_atp` (▽, 9 pts).
  - Notes: Mirrors Fig. 6b's anion-valence ordering: bundling is reversed by `~1.3–1.5 mM` ATP, `~1.8–2.0 mM` ADP, and `~4.5–5.5 mM` AMP.
