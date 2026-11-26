# Figure Digitization Log — Dürre et al. (2018)

Systematic panel-by-panel inventory and calibrated digitization log across all 25 PDF pages of Dürre et al., *Nature Communications* **9**, 1630 (2018) (pp. 1–11 Main Text; pp. 12–25 Supplementary Information).

**Summary**:
- **73 calibrated figure readings** (`figures.yaml` + `figures/*.py`, saved in `readings.json`), totaling **769 digitized data points** across **19 quantitative panels** (10 Main Text + 9 SI), with **0 errors** and **100% visual overlay verification** on contact sheets in `/tmp/figwork_durre2018/`.
- **63 new full-curve PDE reproduction comparisons** (`comparisons.yaml`, 21 predicted model curves × 3 cylinder/cone PDE variants) backed by 6 curve-output wrapper models (`models/*_curves/`) and 12 runs (`runs.yaml`).
- All 52 experimental measurement readings are recorded in `data.yaml` (`origin: measured`) for downstream cross-paper validation; all 100 pre-existing data entries, 13 existing models, refinements, and 21 existing questions are preserved untouched.

---

## Main Text Figures (PDF Pages 1–11)

*Pages 1–2 and 9–11 contain text, Methods, and References only (no figures).*

### Fig. 1 (PDF p. 3)
- **Fig. 1a** — *Left out*: Schematic cartoon of GUV encapsulation components (no quantitative data).
- **Fig. 1b, 1c, 1e** — *Left out*: Confocal fluorescence micrographs of actin inside GUVs (scale bars only; no plotted data).
- **Fig. 1d** — *Read* (**3 readings, 44 pts**): `fig1d_actin` (15 pts, `t = 2..16 min`, `I_norm = -1.1..36.1`), `fig1d_vca_arp_prof` (14 pts, `t = 0.5..16 min`, `31.0..131.5`), `fig1d_cp` (15 pts, `t = 0.5..16 min`, `0.0..109.2`). Rendered at 1200 dpi; calibrated from all inward x/y frame ticks (`0..16 min`, `0..150 a.u.`, linear fit residuals `< 1.0 px`). Darker-series dilation masking prevents cross-series bleed at `t = 2..16 min`. Visually verified.

### Fig. 2 (PDF p. 4)
- **Fig. 2a–g** — *Left out*: Confocal micrographs and 3D fluorescence reconstructions of GUV actin morphologies, VCA/Arp2/3 binding controls, and VCA clustering (no plotted quantitative axes).

### Fig. 3 (PDF p. 5)
- **Fig. 3a–b** — *Left out*: Epifluorescence z-projection time-series micrographs at 40 and 120 nM CP (no plotted quantitative axes).

### Fig. 4 (PDF p. 6)
- **Fig. 4a** — *Read* (**2 readings, 4 pts**): `fig4a_diam_nomyosin` (`60 nM: 4.993 um`, `120 nM: 3.198 um`) and `fig4a_diam_myosin` (`60 nM: 3.172 um`, `120 nM: 3.303 um`), reading horizontal median bars at 1200 dpi (matching caption's `5, 3.2, 3.2 um`). Individual overlapping beeswarm dots are omitted as indistinguishable in dense clusters; y-axis calibrated from ticks `0..10 um` (`res < 1.0 px`). Visually verified.
- **Fig. 4b** — *Read* (**2 readings, 4 pts**): `fig4b_curv_nomyosin` (`60 nM: -0.0992 /um`, `120 nM: -0.1327 /um`) and `fig4b_curv_myosin` (`60 nM: -0.1236 /um`, `120 nM: -0.2203 /um`), matching main text p. 4 (`-0.1, -0.13, -0.12, -0.22 /um`). Calibrated at 1200 dpi from ticks `-0.5..0.0 /um` (`res < 1.0 px`). Visually verified.
- **Fig. 4c** — *Read* (**2 readings, 190 pts**): `fig4c_60nM` (121 light-grey open circle markers, `d = 2.53..9.65 um`, `C = -0.450..-0.016 /um`) and `fig4c_120nM` (69 dark-grey open triangle markers, `d = 1.54..8.11 um`, `C = -0.458..-0.027 /um`), detected via normalized FFT template matching using templates cropped directly from the panel's own legend symbols at 1200 dpi. Calibrated from x-ticks `0..10 um` and y-ticks `-0.5..0.0 /um` (`res < 1.0 px`). Visually verified.
- **Fig. 4d** — *Left out*: Confocal micrographs with zoomed invagination/fission insets (no quantitative axes).

### Fig. 5 (PDF p. 7)
- **Fig. 5a, 5b** — *Left out*: Kinetic pathway schematic (5a) and conceptual 4-regime growth cartoon with unlabelled arbitrary axes (5b).
- **Fig. 5c** — *Read* (**3 readings, 54 pts**): `fig5c_B_40nM` (18 pts, `t = 5..40 s`, peak `444.5 /um^2` at `t = 13 s`, `416.0` at `40 s`), `fig5c_B_80nM` (18 pts, peak `281.1` at `10 s`, `171.2` at `40 s`), `fig5c_B_120nM` (18 pts, peak `198.1` at `10 s`, `97.0` at `40 s`). Calibrated at 1200 dpi from all inward ticks (`0..40 s`, `0..500 /um^2`, `res < 1.0 px`). Compared against all 3 cylinder PDE variants (`cmp_fig5c_*`): printed `0/18` (`not reproduced`, RMSE `142–383`); conserved axisymmetric `not reproduced` (`RMSE 89–133`, peaking `+3.3 s` late and `+29%` high); 1D Cartesian membrane `reproduced` at 40 nM (`18/18`) and `close` at 80/120 nM (`17/18`, RMSE `39–71`). Visually verified.
- **Fig. 5d** — *Read* (**5 readings, 51 pts**): `fig5d_Astar_40nM` (14 pts, peak `78.2` at `4 s`, `4.65` at `40 s`), `fig5d_Astar_80nM` (14 pts, peak `81.6` at `4 s`, `11.77` at `40 s`), `fig5d_Astar_120nM` (17 pts, `t = 1..40 s`, peak `83.5` at `4 s`, `20.06` at `40 s`), plus dashed steady-state guide lines `fig5d_Astar_ss_80nM` (`12.33 /um^2`, 3 pts) and `fig5d_Astar_ss_120nM` (`20.74 /um^2`, 3 pts). Note: y-axis is mislabelled `1/um^2 s` instead of `1/um^2` (`q_fig5d_yaxis_unit_typo`); dashed lines sit `+5–8%` above exact analytical steady states (`11.43, 19.61`, documented in `q_digitize_fig5d_ss_dashed_lines_offset`). Neutral-black tick filtering isolates y-tick `20` cleanly where the dark navy curve crosses the spine (`res < 1.0 px`). Compared in `cmp_fig5d_*`: 1D Cartesian membrane reproduces **5/5** curves (`RMSE 13.2–14.7`, `ss RMSE 0.91–1.13`); axisymmetric fails all 3 transient curves (`RMSE 29.0–34.2`). Visually verified.
- **Fig. 5e** — *Read* (**6 readings, 75 pts**): branching curves `fig5e_branch_40nM` (14 pts, peak `149.4` at `6 s`, `37.43` at `40 s`), `fig5e_branch_80nM` (14 pts, peak `134.6` at `7 s`, `40.73` at `40 s`), `fig5e_branch_120nM` (14 pts, peak `118.4` at `8 s`, `39.03` at `40 s` — directly confirming the non-monotonic `80 > 120 > 40 nM` steady-state ordering in `q_ss_branching_rate_nonmonotonic`); dashed capping curves `fig5e_cap_40nM` (11 pts, peak `56.7`), `fig5e_cap_80nM` (11 pts, peak `80.2`), `fig5e_cap_120nM` (11 pts, peak `85.7`). Calibrated at 1200 dpi (`0..40 s`, `0..160 /(um^2*s)`, `res < 1.0 px`). Compared in `cmp_fig5e_*`: 1D Cartesian membrane reproduces all 3 capping curves (`11/11` each, `RMSE 8.2–21.5`) and scores `close` on all 3 branching curves (`13/14` each, `RMSE 30.7–35.2`); axisymmetric fails all 6 (`RMSE 22.5–76.3`). Visually verified.
- **Fig. 5f** — *Read* (**2 readings, 11 pts**):
  - `fig5f_B_peak` (3 pts at `40, 80, 120 nM`: **`404.2, 254.4, 163.6 /um^2`**), calibrated at 900 dpi directly from the inward ticks on each panel's auto-scaled colorbar (`res <= 1.4 px`). Reveals `120 nM` is `163.6` (not eye-read `~150`); all 3 cone PDE variants score `not reproduced` (`cmp_fig5f_Bpk_*`: peak-1 Cartesian is 8–13% low at `351.3, 225.5, 145.2`, `RMSE 36.4`; axisymmetric is 20–32% high at `486.8, 330.4, 216.3`, `RMSE 71.6`; printed is ~0).
  - `fig5f_80nM_t40_profile` (8 pts at `r = -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, +0.5 um`: **`28.4, 96.0, 131.2, 124.6, 109.5, 99.9, 96.2, 100.3 /um^2`**), calibrated at 1200 dpi on the 3D orthographic front plane (`t = 40 s`) directly from all 5 outward `r`-axis ticks (`-5..+5`, `r_res = 0.20 px`) and 5 `z`-axis ticks (`0..400`, `z_res = 1.29 px`). *Partial exclusion*: `r >= +1.0 um` in 80 nM and the 40/120 nM front edges are left out due to underside 3D wireframe self-occlusion. Reveals a major new discrepancy (`q_digitize_fig5f_80nM_ring_radius_discrepancy`): Fig. 5f crests at `|r| = 2.0 um` and collapses to `28.4 /um^2` at `|r| = 3.0 um`, whereas both conserved cone models crest at `2.91–3.02 um` (`108.4` and `164.9 /um^2` at `r = -3.0 um`, failing `cmp_fig5f_80nM_prof_*`). Visually verified.

### Fig. 6 (PDF p. 8)
- **Fig. 6a** — *Left out*: Confocal micrographs of B16-F1 cells (no quantitative axes).
- **Fig. 6b** — *Read* (**1 reading, 3 pts + whisker range**): `fig6b_microspikes` (`1=CP Dox: 0.00 [0.00, 0.00]`, `2=no CP Dox: 13.00 [0.00, 17.01]`, `3=CP no Dox: 8.00 [0.00, 23.01]`). Calibrated at 1200 dpi (`0..20`, `res < 1.0 px`) with inset micrograph masked and `error_bars: range` (non-outlier whiskers per Fig. 6 caption). Visually verified.
- **Fig. 6c** — *Read* (**1 reading, 3 pts + whisker range**): `fig6c_lamellipodium_thickness` (`1=CP Dox: 0.741 [0.485, 1.042] um`, `2=no CP Dox: 1.480 [1.267, 2.050] um`, `3=CP no Dox: 1.726 [1.030, 2.596] um`), matching text `1.7 -> 0.7 um`. Calibrated at 1200 dpi (`0.5..2.5 um`, `res < 1.0 px`). Visually verified.
- **Fig. 6d** — *Read* (**1 reading, 3 pts + whisker range**): `fig6d_vesicles_per_cell` (`1=CP Dox: 0.291 [0.216, 0.380]`, `2=no CP Dox: 0.126 [0.012, 0.213]`, `3=CP no Dox: 0.119 [0.044, 0.169]`). Calibrated at 1200 dpi (`0.0..0.5`, `res < 1.0 px`, outer whisker-cap wing strip avoiding scatter circle overlap). Quantitatively settles `q_fig6d_vesicle_50pct_vs_2fold` (`q_digitize_fig6d_vesicle_increase_resolved`): `+131%` to `+145%` (`2.31–2.45×`) jump, not `"~50%"`. Visually verified.

---

## Supplementary Information Figures (PDF Pages 12–25)

*Pages 21–25 contain Supplementary Notes 1–2, Supplementary Methods, and Supplementary References (no figures).*

### Supp Fig. 1 (PDF p. 12)
- **Supp Fig. 1** — *Read* (**1 reading, 13 pts**): `suppfig1_pyrene_60nM` (`t = 0..60 min`, `I = -0.135..+0.072 a.u.`, flat baseline confirming bulk inhibition at 60 nM CP). Calibrated at 600 dpi (`0..60 min`, `-1..5 a.u.`, `res < 1.0 px`). Visually verified.

### Supp Fig. 2 (PDF p. 12)
- **Supp Fig. 2 (left/top micrographs)** — *Left out*: Confocal micrographs of GUVs at 0 and 20 nM CP (no quantitative axes).
- **Supp Fig. 2 (bar chart)** — *Read* (**4 readings, 8 pts**): `suppfig2_bulk` (`0 nM: 91.55%`, `20 nM: 21.85%`), `suppfig2_domains` (`4.03%, 36.35%`), `suppfig2_cortex` (`3.98%, 41.51%`), `suppfig2_membrane_total` (`8.39%, 78.14%`). Calibrated at 900 dpi (`0..100%`, `res < 1.0 px`) with Voronoi nearest-color segmentation; matches exact integer counts `N = 24` (`22, 1, 1`) and `N = 96` (`21, 35, 40`) to `< 0.5%`. Visually verified.

### Supp Fig. 3 (PDF p. 13)
- **Supp Fig. 3** — *Left out*: Brightfield GUV morphology micrographs (no quantitative axes).

### Supp Fig. 4 (PDF p. 13)
- **Supp Fig. 4 (3 panels: Domains, Closed Cortex, Bulk Polymerization)** — *Read* (**8 readings, 32 pts** across `c_CP = 40, 60, 120, 180 nM`, `N = 159, 218, 70, 74`):
  - Domains panel: `suppfig4_protrusion` (`24.3, 3.9, 0.0, 1.4%`), `suppfig4_concave_domain` (`32.2, 38.0, 13.7, 0.0%`), `suppfig4_flat_domain` (`33.6, 53.5, 14.6, 7.4%`), `suppfig4_domains_total` (`90.1, 95.5, 28.7, 8.8%`).
  - Closed Cortex panel: `suppfig4_flat_cortex` (`2.9, 0.0, 28.7, 67.7%`), `suppfig4_concave_cortex` (`7.4, 0.5, 31.2, 0.0%`), `suppfig4_cortex_total` (`10.2, 1.0, 59.9, 67.7%`).
  - Bulk Polymerization panel: `suppfig4_bulk` (`0.0, 4.2, 11.5, 24.2%`).
  - Calibrated at 1200 dpi (`0..100%`, `res < 1.0 px`); sums across the 3 panels equal `100.0–100.6%` at every CP concentration. Quantitatively settles `q_closed_cortex_40nM_and_20nM` and `q_audit_reading_quantity_and_180nM_cortex` (`q_digitize_suppfig4_cortex_percentages_resolved`). Visually verified.

### Supp Fig. 5 (PDF p. 14)
- **Supp Fig. 5** — *Left out*: Epifluorescence z-projection time-series micrographs of protrusion formation at 40 nM CP (no quantitative axes).

### Supp Fig. 6 (PDF p. 15)
- **Supp Fig. 6 (inset micrographs)** — *Left out*: Example vesicle micrographs inside the boxplot frame (masked during digitization).
- **Supp Fig. 6 (boxplot)** — *Read* (**1 reading, 3 pts + whisker range**): `suppfig6_domain_count` (`1=protrusion: 15.48 [2.00, 37.99]`, `2=concave: 30.99 [20.00, 37.99]`, `3=flat: 35.04 [8.97, 59.98]`), matching caption's `15.5, 31, 35`. Calibrated at 600 dpi (`0..70`, `res < 1.0 px`) with `error_bars: range`. Visually verified.

### Supp Fig. 7 (PDF p. 16)
- **Supp Fig. 7 (left TIRF micrographs)** — *Left out*: Fluorescence micrograph time-lapse panels (no quantitative axes).
- **Supp Fig. 7 (4 right line-scan panels: 0, 20, 40, 120 nM CP)** — *Read* (**15 readings, 141 pts** across `x = 1..15 um`):
  - `0 nM CP`: `suppfig7_0nM_0min` (8 pts, `x = 1..13`, clipped below y=50 at `x=15`), `1min` (9 pts), `6min` (9 pts), `20min` (9 pts).
  - `20 nM CP`: `suppfig7_20nM_0min`, `1min`, `6min`, `20min` (9 pts each).
  - `40 nM CP`: `suppfig7_40nM_0min` (9 pts, `x = 1..12`, clipped below y=50 at `x=13,15`), `1min`, `3min`, `6min`, `20min` (10 pts each).
  - `120 nM CP`: `suppfig7_120nM_0min` (10 pts), `20min` (10 pts).
  - *Partial exclusions*: Intermediate grey shades that heavily overlap adjacent curves (`3, 9 min` at 0 and 20 nM; `9 min` at 40 nM; `1, 3, 6, 9 min` at 120 nM) are left out as color-indistinguishable. Calibrated at 900 dpi (`0..15 um`, `50..250 a.u.`, `res < 1.0 px`) using darker-curve dilation masking to eliminate anti-aliasing halo pickups. Visually verified.

### Supp Fig. 8 (PDF p. 16)
- **Supp Fig. 8** — *Read* (**5 readings, 58 pts**): `suppfig8_pyrene_20nM` (12 pts, `1..12 min`), `40nM` (9 pts, `1.5..6 min`), `60nM` (8 pts, `2..6 min`), `80nM` (14 pts, `2..14 min`), `180nM` (15 pts, `0..15 min`). *Partial exclusion*: `40` and `60 nM` plateau tails (`t > 6 min`) are left out where they merge identically under the `80 nM` curve. Calibrated at 600 dpi (`0..15 min`, `0..120 A.U.`, `res < 1.0 px`). Visually verified.

### Supp Fig. 9 (PDF p. 17)
- **Supp Fig. 9 (left & right beeswarm plots at 40 nM CP)** — *Read* (**2 readings, 2 pts**): `suppfig9_domain_diameter` (`5.034 um`) and `suppfig9_domain_curvature` (`0.06885 /um`, printed unsigned/positive), matching caption's `~5 um` and `0.07 /um`. Individual overlapping beeswarm dots are omitted as indistinguishable in dense clusters. Calibrated at 1200 dpi (`0..10 um`, `0..0.5 /um`, `res < 1.0 px`). Visually verified.

### Supp Fig. 10 (PDF p. 17)
- **Supp Fig. 10** — *Read* (**4 readings, 8 pts** at `60 nM [N=27]` and `120 nM [N=40]` + 5 uM NMM II): `suppfig10_invagination` (`51.8%, 57.7%`), `suppfig10_fission` (`3.9%, 27.4%`), `suppfig10_deformation` (`29.4%, 7.7%`), `suppfig10_no_deformation` (`14.8%, 7.6%`), matching exact integer counts (`14, 1, 8, 4 / 27` and `23, 11, 3, 3 / 40`) to `< 0.5%`. Calibrated at 900 dpi (`0..100%`, `res < 1.0 px`). Visually verified.

### Supp Fig. 11 (PDF p. 18)
- **Supp Fig. 11** — *Read* (**3 readings, 57 pts** across `t = 6..100 s`): `suppfig11_B_40nM` (19 pts, `t = 100 s: 468.1 /um^2`), `suppfig11_B_80nM` (19 pts, `176.0 /um^2`), `suppfig11_B_120nM` (19 pts, `98.6 /um^2`). Calibrated at 900 dpi (`0..100 s`, `0..500 /um^2`, `res < 1.0 px`). Compared in `cmp_sfig11_*`: 1D Cartesian membrane reproduces **3/3** curves (`19/19` each, `RMSE 26.6–52.8`); axisymmetric fails transients (`RMSE 66.1–93.7`) while matching `t = 100 s` steady states (`471.5, 176.8, 98.7`) to `< 0.8%`; printed fails completely (`RMSE 131–342`). Visually verified.

### Supp Fig. 12 (PDF p. 18)
- **Supp Fig. 12** — *Left out*: Western blot gel image of CP alpha1/beta2 and tubulin in B16-F1 cells (no plotted quantitative axes; densitometry `2.1-fold` already recorded in `data.yaml`).

### Supp Fig. 13 (PDF p. 19)
- **Supp Fig. 13** — *Left out*: Confocal immunofluorescence colocalization micrographs of CP, F-actin, and LAMP1/EEA1/Rab7 in B16-F1 cells (no quantitative axes).

### Supp Fig. 14 (PDF p. 20)
- **Supp Fig. 14** — *Read* (**2 readings, 6 pts** across `40, 80, 120 nM CP`):
  - `suppfig14_Vtot_max` (**`2131, 1928, 1575 /um^2`**) and `suppfig14_Vtot_min` (**`819.5, 716.3, 651.1 /um^2`**), calibrated at 900 dpi from inward colorbar ticks via saturation-drop peak detection (`res <= 2.3 px`).
  - *Partial exclusion*: The 3D front-plane VCA surface edge (`t = 40 s`) is left out because the back-wall surface protrudes beneath the arched foreground mesh across `r >= 0`.
  - Yields exact `t = 40 s` max/min VCA enrichment ratios of **`2.60×, 2.69×, 2.42×`**, quantitatively confirming `q_vca_4fold_time_window`. Compared in `cmp_sfig14_*`: 1D Cartesian membrane reproduces **both** `V_tot_max` (`3/3`, `RMSE 99.7`) and `V_tot_min` (`3/3`, `RMSE 40.1`) within input rounding; axisymmetric overpredicts both (`0/3`, `RMSE 251.2` and `116.4`); printed fails (`0/3`, `RMSE 820.8` and `350.2`). Visually verified.
