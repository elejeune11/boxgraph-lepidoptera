# Figure Digitization Log: `bashirzadeh2021`

Complete calibrated figure inventory and digitization log for Bashirzadeh, Redford, Lorpaiboon, Groaz, Moghimianavval, Litschel, Schwille, Hocky, Dinner & Liu, *Commun. Biol.* **4**, 1136 (2021), covering all **31 PDF pages** (Main Text pp. 1–11; Supplemental Information pp. 12–31 / SI-1–SI-20).

All **25 quantitative figure panels** across 8 pages (pp. 3, 4, 5, 6, 17, 21, 23, 25) were digitized at **900–1200 dpi** in `mode: code` using self-calibrating readers in `figures/` (shared fast vectorized engine [`figures/_calib.py`](file:///home/antigravity/boxgraph-public/corpora/actin_bundling/records/bashirzadeh2021/figures/_calib.py), 11–189 ms per reading directly on the passed crop array), yielding **60 calibrated figure readings** (**495 data points**, **0 missing points**, **21 readings with digitized error bars**) and **22 new quantitative comparisons** (46 total comparisons). Every reading was run through `boxgraph.engine.Builder` / `boxgraph figure check`, saved into `readings.json`, and visually verified across all 60 overlays on compact 3×2 contact sheets (including high-contrast error-bar cap verification sheets) in `/tmp/figwork_bashirzadeh2021/`.

---

## Main Article Figures (PDF pp. 3–7)

### Main Fig. 1 (p. 3) — GUV size- and α-actinin concentration-dependent organization
- **Fig. 1a** *(left out)*: Representative 3D-reconstructed fluorescence and skeletonized confocal micrographs of α-actinin-crosslinked actin networks in small, medium, and large GUVs. Qualitative microscopy images with no numeric plot axes.
- **Fig. 1b** *(read — 2 readings, 6 bars total)*:
  - **`fig1b_ring_prob`**, **`fig1b_ring_plus_ringnet_prob`** (`fig1b_*_fig`, `measured`, 900 dpi, crop `[1040, 2240, 2310, 3530]`): Stacked cumulative probability bars (Ring, Ring/Network, Network) across `7-12`, `12-16`, `>16 um` GUV bins calibrated on 6 outward y-ticks (`0.0..1.0`, residual `0.15 px`).
  - Readings: `P(Ring) = [0.962, 0.506, 0.108]`; `P(Ring + Ring/Network) = [0.983, 0.875, 0.590]` (implying `P(Network) = [0.017, 0.125, 0.410]`).
  - **Comparison**: `cmp_fig1b_ring_vs_nguv_weighted_fig1c` compares `fig1b_ring_prob_fig` against the `NGUVs`-weighted average of Fig. 1c (`[0.956, 0.465, 0.039]`): large (`>16 um`) GUVs diverge by `> 2.5×` (`0.108` vs `0.039`, and also exceed the unweighted average `0.079`), logged in `q_digitize_fig1b_fig2b_pooled_vs_mean`.
- **Fig. 1c** *(read — 3 readings, 9 bars with ±SEM error bars)*:
  - **`fig1c_ring_pct_0p1`**, **`fig1c_ring_pct_0p2`**, **`fig1c_ring_pct_0p3`** (`fig1c_ring_pct_*_fig`, `measured`, 900 dpi, crop `[2370, 2230, 3540, 3540]`, `error_bars: sem`): Hollow colored bars (`α-actinin/actin = 0.1, 0.2, 0.3 M/M`) across `7-12, 12-16, >16 um` calibrated on 11 outward y-ticks (`0..100%`, residual `0.53 px`).
  - Readings (mean `[bar_lo, bar_hi]`):
    - `0.1`: `99.8%` (no bar), `55.7% [50.0, 61.5]`, `18.7% [11.6, 25.9]`
    - `0.2`: `90.5% [81.0, 100.0]`, `54.5% [40.9, 68.3]`, `3.6% [-0.1, 7.3]`
    - `0.3`: `88.8% [81.0, 96.5]`, `41.3% [32.7, 49.9]`, `1.3%` (occluded by replicate diamonds)
- **Fig. 1d** *(read — 6 readings, 18 stacked bar heights total)*:
  - **`fig1d_network_pct_0p1/0p2/0p3`** (bottom `Network` layer) and **`fig1d_total_network_pct_0p1/0p2/0p3`** (stacked `Network + Ring/Network` top) (`fig1d_*_fig`, `measured`, 900 dpi, crop `[3540, 2230, 4740, 3540]`): Calibrated on 7 outward y-ticks (`0..120%`, residual `0.32 px`). No error bars are drawn (despite the caption's typo `"Error bars in (d)"`, already logged in `q_figure_labels_and_stacked_diamonds_gt100`).
  - Readings (`Network` / `Network + Ring/Network`):
    - `0.1`: `0.0% / 0.0%`, `0.0% / 43.6%`, `19.2% / 80.6%`
    - `0.2`: `0.0% / 9.5%`, `3.9% / 45.4%`, `32.1% / 97.2%`
    - `0.3`: `4.8% / 11.1%`, `23.6% / 58.8%`, `61.7% / 98.9%`
  - **Consistency checks**: `cmp_fig1cd_partition_0p1/0p2/0p3` verify `Ring% + (Network + Ring/Network)% = 100%` across all 9 `(ratio, size)` bins (`99.3%–100.8%`, every bin within `0.8%` of `100%`).
- **Fig. 1e** *(left out)*: 3D-reconstructed and skeletonized confocal micrographs of peripheral asters in large GUVs. No numeric axes.
- **Fig. 1f** *(read — 3 readings, 6 bars with one-sided upward error whiskers)*:
  - **`fig1f_periph_aster_pct_0p1/0p2/0p3`** (`fig1f_periph_aster_pct_*_fig`, `measured`, 900 dpi, crop `[1050, 4140, 2340, 5550]`, `error_bars: sem` per caption): Calibrated on 6 outward y-ticks (`0..100%`, residual `0.40 px`) across `7-16` and `>16 um`.
  - Readings (bar top / upper whisker endpoint): `0.1`: `8.3% (13.3%)`, `53.3% (67.1%)`; `0.2`: `12.9% (26.7%)`, `53.0% (60.0%)`; `0.3`: `17.3% (35.6%)`, `64.9% (88.3%)`.
  - **Discrepancy discovered (`q_digitize_fig1f_fig2f_onesided_sd_vs_sem`)**: Although the caption states `"Error bars indicate standard error of the mean; n = 3 experiments"`, all stems are uncapped one-sided upward whiskers terminating at the highest of the 3 replicate diamonds (`~mean + 1 SD`, e.g., `64.9% -> 88.3%` vs `mean + SEM = 78.4%`), not two-sided T-capped `±SEM` bars as in Fig. 1c.
- **Fig. 1g** *(read — 1 reading, 3 stacked bars)*:
  - **`fig1g_periph_bundles_pct`** (`fig1g_periph_bundles_pct_fig`, `measured`, 900 dpi, crop `[2370, 4230, 3540, 5550]`): Percentage of peripheral asters with peripheral bundles (vs lumenal bundles) across `α-actinin/actin = [0.1, 0.2, 0.3]`, calibrated on 6 outward y-ticks (`0..100%`, residual `0.25 px`). Readings: **`43.0%`, `58.5%`, `67.9%`** (implying lumenal bundles `57.0%, 41.5%, 32.1%`).
- **Fig. 1h** *(left out)*: Schematic summary diagram of α-actinin-actin structures as a function of GUV size and crosslinker concentration. No numeric data.

---

### Main Fig. 2 (p. 4) — α-actinin + fascin (0.1 M/M) actin structures in GUVs
- **Fig. 2a** *(left out)*: 3D-reconstructed and skeletonized confocal micrographs of rings, peripheral asters, and central asters. No numeric axes.
- **Fig. 2b** *(read — 2 readings, 6 bars total)*:
  - **`fig2b_ring_prob`**, **`fig2b_ring_plus_ringnet_prob`** (`fig2b_*_fig`, `measured`, 900 dpi, crop `[360, 3600, 1850, 5020]`): Calibrated on 6 outward y-ticks (`0.0..1.0`, residual `0.20 px`) across `7-11` (printed label for `7-12`), `12-16`, `>16 um`.
  - Readings: `P(Ring) = [0.779, 0.120, 0.019]`; `P(Ring + Ring/Network) = [0.983, 0.927, 0.634]` (implying `P(Network) = [0.017, 0.073, 0.366]`).
  - **Comparison**: `cmp_fig2b_ring_vs_nguv_weighted_fig2d` reveals that in `12-16 um` GUVs, cumulative Fig. 2b `P(Ring) = 0.120` falls strictly below EVERY individual Fig. 2d ratio mean (`0.131, 0.198, 0.203`; `NGUVs`-weighted `0.172`, unweighted `0.178`), documenting Simpson-type per-experiment weighting or pooling inconsistency (`q_digitize_fig1b_fig2b_pooled_vs_mean`).
- **Fig. 2c** *(left out)*: Conceptual schematic of GUV-size-dependent α-actinin+fascin architectures. No numeric data.
- **Fig. 2d** *(read — 3 readings, 9 bars with ±SEM error bars)*:
  - **`fig2d_ring_pct_0p1/0p2/0p3`** (`fig2d_ring_pct_*_fig`, `measured`, 900 dpi, crop `[360, 5020, 1850, 6500]`, `error_bars: sem`): Calibrated on 11 outward y-ticks (`0..100%`, residual `0.82 px`).
  - Readings (mean `[bar_lo, bar_hi]`):
    - `0.1`: `62.4% [53.8, 71.0]`, `13.1%` (bar top beneath `15.5%` diamond; SEM hidden by diamonds), `1.2%` (near-zero diamond at `0%`)
    - `0.2`: `79.5% [74.8, 84.3]`, `19.8% [6.1, 33.5]`, `5.6% [2.7, 8.5]`
    - `0.3`: `78.4% [75.0, 81.8]`, `20.3% [10.4, 30.3]`, `1.2%` (near-zero diamond at `0%`)
- **Fig. 2e** *(read — 6 readings, 18 stacked bar heights total)*:
  - **`fig2e_network_pct_0p1/0p2/0p3`** (bottom `Network`) and **`fig2e_total_network_pct_0p1/0p2/0p3`** (stacked `Network + Ring/Network`) (`fig2e_*_fig`, `measured`, 900 dpi, crop `[1820, 5020, 3260, 6500]`): Calibrated on 7 outward y-ticks (`0..120%`, residual `0.46 px`).
  - Readings (`Network` / `Network + Ring/Network`):
    - `0.1`: `5.0% / 38.6%`, `22.7% / 86.7%`, `52.1% / 100.1%`
    - `0.2`: `0.0% / 20.0%`, `0.0% / 78.8%`, `18.3% / 95.9%`
    - `0.3`: `0.0% / 20.4%`, `6.1% / 79.3%`, `37.9% / 99.8%`
  - **Consistency checks**: `cmp_fig2de_partition_0p1/0p2/0p3` confirm all 9 `(ratio, size)` bins sum to `98.7%–101.5%` (~100%).
- **Fig. 2f** *(read — 3 readings, 6 bars with one-sided upward error whiskers)*:
  - **`fig2f_aster_pct_0p1/0p2/0p3`** (`fig2f_aster_pct_*_fig`, `measured`, 900 dpi, crop `[3220, 5020, 4700, 6500]`, `error_bars: sem` per caption): Calibrated on 6 outward y-ticks (`0..100%`, residual `0.35 px`).
  - Readings (bar top / upper whisker endpoint): `0.1`: `24.1% (25.9%)`, `65.8% (78.8%)`; `0.2`: `15.5% (20.3%)`, `61.3% (66.8%)`; `0.3`: `19.9% (20.4%)`, `61.3% (77.1%)`. Shows the same one-sided uncapped `~mean + 1 SD` whisker discrepancy as Fig. 1f (`q_digitize_fig1f_fig2f_onesided_sd_vs_sem`).

---

### Main Fig. 3 (p. 5) — GUV-size-dependent F-actin cluster localization
- **Fig. 3a** *(left out)*: Representative 2D and 3D-reconstructed confocal fluorescence micrographs of central asters at `α-actinin/actin = 0.3` and `0.1` (with `fascin/actin = 0.1`). No numeric axes.
- **Fig. 3b** *(read — 2 readings, 37 points total)*:
  - **`fig3b_actin_intensity_0p3`** (black, 18 pts `x = 0.05..0.90`), **`fig3b_actin_intensity_0p1`** (red, 19 pts `x = 0.05..0.98`) (`fig3b_actin_intensity_*_fig`, `measured`, 900 dpi, crop `[360, 2650, 2100, 4120]`): Calibrated on 6 major outward x-ticks (`0.0..1.0`, residual `0.31 px`) and 6 outward y-ticks (`0..200 a.u.`, residual `0.35 px`). Captures the broad bimodal high-intensity central cluster plateau (`~185 a.u.` over `x = 0.32..0.58`) at `0.3` vs the narrow central spike (`~150 a.u.` at `x = 0.50`) at `0.1`.
- **Fig. 3c** *(read — 2 readings, 4 bars total)*:
  - **`fig3c_cluster_center_pct`**, **`fig3c_cluster_periphery_pct`** (`fig3c_cluster_*_fig`, `measured`, 900 dpi, crop `[2080, 2650, 3450, 4120]`): Calibrated on 6 outward y-ticks (`0..100%`, residual `0.28 px`) and 2 outward x-ticks (`α-Actinin`, `α-Actin+fascin` [sic]).
  - Readings (`Center` / `Periphery` / `Sum`): α-actinin alone (`>20 um`): **`3.7%` / `62.4%`** (sum `66.1%`, `33.9%` unclustered); α-actinin + fascin (`>20 um`): **`29.8%` / `47.1%`** (sum `76.9%`, `23.1%` unclustered).
- **Fig. 3d** *(read — 2 readings, 4 bars total)*:
  - **`fig3d_cluster_center_pct`**, **`fig3d_cluster_periphery_pct`** (`fig3d_cluster_*_fig`, `measured`, 900 dpi, crop `[3420, 2650, 4740, 4120]`): Calibrated on 6 outward y-ticks (`0..100%`, residual `0.30 px`) and 2 outward x-ticks (`16-24`, `>24 um`).
  - Readings (`Center` / `Periphery` / `Sum`): `16-24 um`: **`18.0%` / `55.6%`** (sum `73.6%`); `>24 um`: **`57.7%` / `26.9%`** (sum `84.6%`).

---

### Main Fig. 4 (p. 6; caption p. 7) — Crosslinker sorting, bundle persistence length $L_p$, and AFINES simulations
- **Fig. 4a, 4b** *(left out)*: Representative confocal micrographs of TMR α-actinin, AF647 fascin, and ATTO-488 actin in encapsulated central asters. No numeric axes.
- **Fig. 4c** *(read — 2 readings, 37 points total)*:
  - **`fig4c_fascin_intensity`** (magenta, 19 pts `x = 1..47 um`), **`fig4c_alpha_actinin_intensity`** (blue, 18 pts `x = 3..45 um`) (`fig4c_*_fig`, `measured`, 1200 dpi, crop `[4700, 1280, 6950, 2900]`): Calibrated on 4 inward y-ticks (`0.5..2.0`, residual `0.34 px`) and 3 x-axis references (`0, 25, 50 um`). Shows α-actinin strongly enriched at the central cluster (`~1.83` at `x = 24 um`) vs broader fascin distribution (`~1.35`).
- **Fig. 4d** *(read — 2 readings, 36 points total)*:
  - **`fig4d_fascin_angular_intensity`** (magenta, left axis `1400..1900 a.u.`, 18 pts `20..360 deg`), **`fig4d_alpha_actinin_angular_intensity`** (blue, right axis `500..1000 a.u.`, 18 pts `20..360 deg`) (`fig4d_*_fig`, `measured`, 1200 dpi, crop `[7050, 1250, 9650, 2900]`): Dual-axis calibration on 6 inward y-ticks (residual `0.42 px`) and 3 inward x-ticks (`0, 200, 400 deg`).
- **Fig. 4e** *(read — 2 readings, 5 bars with symmetric ±SEM error bars)*:
  - **`fig4e_lp_no_fascin`** (`[0:1:10]`, `[0:3:10]`), **`fig4e_lp_plus_fascin`** (`[1:1:10]`, `[1:2:10]`, `[1:3:10]`) (`fig4e_lp_*_fig`, `measured`, 1200 dpi, crop `[1100, 3150, 4250, 5600]`, `error_bars: sem`): Calibrated on 6 intensity-weighted outward y-ticks (`0..250 um`, residual `0.95 px`) and 5 outward x-ticks.
  - Readings (bar mean `±` plotted SEM half-width vs caption printed mean `±` value):
    - `[0:1:10]` (`n = 22`): **`33.60 ± 7.72 um`** (`[25.88, 41.33]`) vs caption `33.3 ± 10 um` (**−23% narrower SEM**)
    - `[0:3:10]` (`n = 14`): **`39.40 ± 8.69 um`** (`[30.71, 48.09]`) vs caption `39.3 ± 9 um`
    - `[1:1:10]` (`n = 17`): **`54.36 ± 10.38 um`** (`[43.74, 64.50]`) vs caption `54.2 ± 11.4 um` (**−9% narrower SEM**)
    - `[1:2:10]` (`n = 14`): **`78.01 ± 19.07 um`** (`[58.71, 96.84]`) vs caption `77.8 ± 19 um`
    - `[1:3:10]` (`n = 26`): **`99.74 ± 14.00 um`** (`[85.74, 113.74]`) vs caption `99.7 ± 15 um` (**−7% narrower SEM**)
  - **Key quantitative discovery — resolves `q_welch_pvalues_fig4e` (`q_digitize_fig4e_errorbars_resolve_welch`)**:
    - `cmp_fig4e_bar_*_vs_caption` (5/5 **`reproduced`**) confirm all 5 bar heights match the caption means within `0.04–0.30 um`.
    - `cmp_fig4e_p_0110_vs_1110_plotted_sems`, `cmp_fig4e_p_1110_vs_1310_plotted_sems`, and `cmp_fig4e_p_0310_vs_1310_plotted_sems` (3/3 **`reproduced`**!) show that running Welch's two-sample t-test with the **plotted Fig. 4e SEMs** reproduces ALL THREE printed p-values to `~0.00005`: **`0.11625`** vs `0.1162` (caption SEMs gave `0.1770`), **`0.01258`** vs `0.0125` (caption SEMs gave `0.0203`), and **`0.00076`** vs `0.0008` (caption SEMs gave `0.00140`).
- **Fig. 4f** *(left out)*: Final-frame 2D AFINES simulation snapshot (actin red, α-actinin cyan, fascin black). Qualitative spatial snapshot.
- **Fig. 4g** *(left out)*: Four 20×20 2D spatial probability density heatmaps of α-actinin and fascin (single vs combined) with colorbar spans `0..3×10^-3`, `0..5×10^-3`, `0..1.25×10^-2 um^-2`. Qualitative spatial heatmaps without extractable scalar summary curves (normalization already checked in `run_hist_pdf_uniform_disk`).
- **Fig. 4h** *(read — 4 readings, 40 points with symmetric ±SD error bars)*:
  - **`fig4h_single_alpha_actinin`**, **`fig4h_single_fascin`** (top panel, crop `[4600, 4430, 7400, 6220]`) and **`fig4h_combined_alpha_actinin`**, **`fig4h_combined_fascin`** (bottom panel, crop `[4600, 6220, 7400, 8150]`) (`fig4h_*_fig`, `run`, 1200 dpi, `error_bars: sd`): Box-frame calibrated on top spine (`1.000`) + 6 inward y-ticks (`0.970..0.995`, residual `0.35 px`) and 4 interior bottom x-ticks (`20, 40, 60, 80 s`, residual `0.40 px`), with all 10 evenly spaced time stations (`t = 3.5–94.4 s`) isolated directly from the red/black vertical SD error bars.
  - Readings:
    - Single α-actinin: rises from `0.9802 ± 0.0039` (`3.5 s`) to `0.9945–0.9952` (`23.7–94.4 s`); Single fascin: flat `0.9950–0.9961`.
    - Combined α-actinin: rises from `0.9770 ± 0.0059` (`3.6 s`) to `0.9852–0.9885` (`±0.001–0.003`); Combined fascin: flat `0.9947–0.9961`.
  - **Model comparisons**: `cmp_fig4h_single_alpha_vs_eq5_bound`, `cmp_fig4h_single_fascin_vs_eq5_bound`, `cmp_fig4h_combined_alpha_vs_eq5_bound`, and `cmp_fig4h_combined_fascin_vs_eq5_bound` test all **40/40 digitized points** against printed Eq. 5's unconditional Cauchy–Schwarz ceiling $U_4 \le 2/3 = 0.6667$ (`binder_cumulant_eq5_bound_vs_time`); **40/40 points** are **`not reproduced`**, exceeding $2/3$ by `+0.310` to `+0.329` (`q_digitize_fig4h_trajectories_vs_eq5`).
- **Fig. 4i** *(left out)*: Conceptual cartoon schematic of crosslinker sorting during central aster assembly. No numeric data.

---

## Supplemental Information Figures (PDF pp. 15–28 / SI-4–SI-17)

- **Supplementary Fig. 1 (p. 15, SI-4)** *(left out)*: cDICE schematic (a) and confocal micrographs of actin-only and α-actinin-actin GUVs (b, c). No quantitative plots.
- **Supplementary Fig. 2 (p. 16, SI-5)** *(left out)*: Image-processing/SOAX skeletonization workflow chart (a) and 3D/skeletonized GUV micrographs (b). No quantitative plots.
- **Supplementary Fig. 3 (p. 17, SI-6)**:
  - **Supp. Fig. 3a, 3b** *(left out)*: 3D/skeletonized actin ring micrographs (a) and geometric schematic defining $d$, $d_{\text{avg}}$, $s$, $p$ (b).
  - **Supp. Fig. 3c** *(read — 2 readings, 32 points total)*: **`supp_fig3c_ring_dev_small`** (blue, 17 pts `s/p = 0.04..0.84`), **`supp_fig3c_ring_dev_medium`** (red, 15 pts `s/p = 0.04..0.78`) (`supp_fig3c_ring_dev_*_fig`, `measured`, 1200 dpi, crop `[5500, 7250, 9400, 11100]`): Calibrated on 3 inward x-ticks (`0, 0.5, 1`) and 3 y-axis ticks (`0.5, 1.0, 1.5`). Quantifies periodic radial deviations $d/d_{\text{avg}} \in [0.91, 1.08]$ along skeletonized small and medium rings.
- **Supplementary Fig. 4 (p. 18, SI-7)** *(left out)*: 3D confocal micrographs (a, c) and peripheral-aster schematic (b). No numeric axes.
- **Supplementary Fig. 5 (p. 19, SI-8)** *(left out)*: Confocal micrographs of fascin-bundled actin and sfGFP-fascin in GUVs (a, b). No numeric axes.
- **Supplementary Fig. 6 (p. 20, SI-9)** *(left out)*: Time-lapse confocal fluorescence micrographs of central aster assembly (`0–30 min`). No numeric axes.
- **Supplementary Fig. 7 (pp. 21–22, SI-10–SI-11)**:
  - **Supp. Fig. 7a, 7c** *(left out)*: Confocal micrographs of bulk fascin-actin bundles (a) and encapsulated central asters (c).
  - **Supp. Fig. 7b** *(read — 1 reading, 22 points)*: **`supp_fig7b_bundle_profile`** (`supp_fig7b_bundle_profile_fig`, `measured`, 1200 dpi, crop `[6150, 1950, 9400, 5200]`): Single-bundle fluorescence intensity profile across `0..10 um` (inward y-ticks `100, 150, 200 a.u.`, x-ticks `0, 5, 10 um`), cleanly masking the red dashed `Δ` annotation and capturing the sharp bundle peak (`~192 a.u.` at `x = 4.6 um`, baseline `~100 a.u.`).
  - **Supp. Fig. 7d** *(read — 2 readings, 4 bars with ±SD error bars)*: **`supp_fig7d_bundle_width_bulk`** (`fascin/actin = 0.05, 0.1, 0.5`), **`supp_fig7d_bundle_width_aster`** (`0.1`, red bar) (`supp_fig7d_bundle_width_*_fig`, `measured`, 1200 dpi, crop `[5500, 6250, 9450, 10150]`, `error_bars: sd`): Calibrated on 7 outward y-ticks (`0.0..3.0 um`, residual `0.38 px`) and 4 outward x-ticks.
    - Readings: Bulk `0.05`: **`0.841 ± 0.181 um`**; Bulk `0.1`: **`1.278 ± 0.247 um`**; Bulk `0.5`: **`1.864 ± 0.733 um`**; Encapsulated central aster (`0.1`, `α-actinin/actin = 0.3`): **`1.173 ± 0.277 um`**.
    - **Comparison**: `cmp_supp_fig7d_bulk_0p1_vs_aster` (**`reproduced`**) confirms bulk `0.1` (`1.278 um`) lies cleanly inside the central aster SD band (`[0.896, 1.450] um`), quantitatively corroborating Main Text p. 5.
- **Supplementary Fig. 8 (pp. 23–24, SI-12–SI-13)**:
  - **Supp. Fig. 8a, 8d** *(left out)*: Confocal/skeletonized fluorescence micrographs of star-like networks and actin rings.
  - **Supp. Fig. 8b** *(read — 2 readings, 36 points total)*: **`supp_fig8b_actin_angular`** (green, 18 pts), **`supp_fig8b_alpha_actinin_angular`** (magenta, 18 pts) (`supp_fig8b_*_fig`, `measured`, 900 dpi, crop `[1150, 3950, 4150, 6000]`): Normalized angular intensity (`0.5..2.0`, `0..360 deg`) around a star-like network outside the central cluster.
  - **Supp. Fig. 8c** *(read — 2 readings, 42 points total)*: **`supp_fig8c_actin_ring_profile`** (green, 21 pts), **`supp_fig8c_alpha_actinin_ring_profile`** (magenta, 21 pts) (`supp_fig8c_*_fig`, `measured`, 900 dpi, crop `[4250, 3950, 7250, 6000]`): Normalized line scan (`0..15 um`, `0.0..2.0`) across an actin ring showing dual wall peaks near `x ≈ 1.0` and `11.7 um`.
  - **Supp. Fig. 8e** *(read — 3 readings, 57 points total)*: **`supp_fig8e_actin_ring_profile`** (green), **`supp_fig8e_alpha_actinin_ring_profile`** (magenta), **`supp_fig8e_fascin_ring_profile`** (blue) (19 pts each, `x = 0.2..3.8 um`) (`supp_fig8e_*_fig`, `measured`, 1200 dpi, crop `[5500, 8250, 9150, 12350]`): Normalized ring cross-section (`0.8..1.4`, `0..4 um`) showing colocalized actin (`~1.48`), α-actinin (`~1.23`), and fascin (`~1.07`) peaks near `x = 2.2–2.4 um`.
- **Supplementary Fig. 9 (p. 25, SI-14)**:
  - **Supp. Fig. 9a** *(left out)*: 2D and 3D confocal micrographs of central asters at `fascin = 0.5 uM` (`0.1`) and `1.5 uM` (`0.3`).
  - **Supp. Fig. 9b** *(read — 2 readings, 39 points total)*: **`supp_fig9b_actin_intensity_0p3`** (black, 19 pts `x = 1..30 um`, sampling unoccluded stations around the `x = 18` overlap), **`supp_fig9b_actin_intensity_0p1`** (red, 20 pts `x = 1..30 um`) (`supp_fig9b_*_fig`, `measured`, 1200 dpi, crop `[1150, 6450, 4950, 10850]`): Calibrated on 6 outward y-ticks (`0..200 a.u.`) and 9 outward x-ticks (`-5..35 um`), showing comparable bimodal central cluster peak intensities (`~190 a.u.`).
  - **Supp. Fig. 9c** *(read — 1 reading, 2 bars with ±SEM error bars)*: **`supp_fig9c_bundle_lp`** (`supp_fig9c_bundle_lp_fig`, `measured`, 1200 dpi, crop `[5300, 6450, 9100, 10850]`, `error_bars: sem`): Calibrated on 7 outward y-ticks (`0..120 um`, residual `0.42 px`), with bar centers located directly from the 4 vertical bar side walls (bypassing the author's off-center x-tick under `"0.3"`).
    - Readings: `fascin/actin = 0.1` (`n = 26`, identical `[1:3:10]` condition as Fig. 4e Bar 5): **`99.58 ± 15.09 um`** (`[84.60, 114.78]`); `fascin/actin = 0.3` (`n = 16`): **`97.66 ± 21.79 um`** (`[75.76, 119.34]`).
    - **Cross-figure finding (`cmp_supp_fig9c_bar1_vs_lp_1_3_10`, `q_digitize_fig4e_errorbars_resolve_welch`)**: Supp. Fig. 9c Bar 1 matches the Fig. 4e caption's `99.7 ± 15 um` BOTH in mean (`99.58 um`, **`reproduced`**) AND in SEM (`±15.09 um`), directly contrasting with Fig. 4e Bar 5's narrower plotted SEM (`±14.00 um`).
- **Supplementary Fig. 10 (p. 26, SI-15)** *(left out)*: Four 2D AFINES simulation snapshots (initial condition + 100 s fascin-only, α-actinin-only, and combined states). Qualitative snapshots with no numeric axes.
- **Supplementary Fig. 11 (p. 27, SI-16)** *(left out)*: Time-lapse 2D count-histogram heatmaps (`t = 0..100 s`) of bound α-actinin (single vs combined) without numeric color-scale summary curves.
- **Supplementary Fig. 12 (p. 28, SI-17)** *(left out)*: Two Cytosim 50 s simulation snapshots with and without the fork feature (color-legend vs caption swap already logged in `q_supp_fig12_legend_caption_swap`). No numeric axes.
