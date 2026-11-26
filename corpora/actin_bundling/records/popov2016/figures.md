# Figure Digitization Log: `popov2016`

Complete calibrated figure inventory and digitization log for Popov, Komianos & Papoian, *PLoS Comput. Biol.* **12**(4), e1004877 (2016), covering all 62 PDF pages (Main Text pp. 1–35; Supporting Information S1–S4 Texts and S1–S2 Tables pp. 36–62).

All 12 quantitative figure panels across 10 pages were digitized at **900 dpi** in `mode: code` (using self-calibrating readers in `figures/` operating directly on the passed crop array), yielding **39 calibrated figure readings** (345 data points) and **23 new comparisons**. Every reading was verified with `boxgraph figure check` and visually inspected on overlay contact sheets in `figwork/`.

---

## Main Article Figures (PDF pp. 7–26)

### Main Fig. 1 (p. 7) — Cytoskeletal network schematic
- **Fig. 1** *(left out)*: 3D schematic cartoon of a MEDYAN cytoskeletal network (Arp2/3, Formin, Ena/VASP, Capping Protein, Actin, NM II Motor, Filament Bundling on a 100 nm grid). No quantitative data or numeric axes.

### Main Fig. 2 (p. 9) — Cylindrical polymer segment schematic
- **Fig. 2** *(left out)*: Geometric schematic of coarse-grained cylindrical polymer segments illustrating $\sigma_0$, $l_i$, $\delta l \ll l_0$, and $\theta_{i,i+1}$. No quantitative data.

### Main Fig. 3 (p. 11) — Mechanical interaction schematics
- **Fig. 3A–C** *(left out)*: Three geometric schematics showing excluded-volume coordinate parameterization (A), branching junction geometry $\theta_{i,j}, d_{i,j}$ (B), and motor binding sites $x_i^m, x_j^m$ (C). No quantitative data.

### Main Fig. 4 (p. 15) — Simulation workflow diagram
- **Fig. 4** *(left out)*: Four-step conceptual flow diagram of the MEDYAN mechanochemical simulation loop. No quantitative data.

### Main Fig. 5 (p. 18) — Single trajectory 3D snapshot
- **Fig. 5** *(left out)*: Single trajectory 3D rendering of a $1 \times 1 \times 1\ \mu\text{m}^3$ actomyosin network at $R_{\alpha:a} = 0.1$, $R_{m:a} = 0.01$ after $2000\ \text{s}$. Qualitative structural snapshot with no numeric axes.

### Main Fig. 6 (p. 19) — Heat map of network $R_g$ ($\mu$m) after 2000 s
- **Fig. 6** *(read — 3 readings, 18 cells total)*:
  - **`fig6_rg_rma0p005`**, **`fig6_rg_rma0p01`**, **`fig6_rg_rma0p02`** (`fig6_rg_rma*_fig`, `run`, `figures/fig6_rg_rma*.py` via `_heatmaps.py`, 900 dpi, crop `[2400, 820, 6900, 4340]`): All three rows ($R_{m:a} = 0.005, 0.01, 0.02$) across all six columns ($R_{\alpha:a} = 0.01, 0.02, 0.05, 0.1, 0.2, 0.5$), decoded directly from the image's own vertical colorbar calibrated on its five printed labels `1.8, 1.6, 1.4, 1.2, 1.0` (spacing `720.0 px` per `0.2 um`, calibration residual `0.25 px`, exact span `[1.000, 1.800] um`).
  - Readings:
    - $R_{m:a} = 0.005$: `1.797, 1.715, 1.532, 1.423, 1.449, 1.382 um`
    - $R_{m:a} = 0.01$: `1.645, 1.535, 1.355, 1.283, 1.257, 1.197 um`
    - $R_{m:a} = 0.02$: `1.464, 1.361, 1.100, 1.046, 1.004, 1.005 um`
  - **Model comparisons**: `cmp_fig6_rma0p005_vs_1um_bound`, `cmp_fig6_rma0p01_vs_1um_bound`, `cmp_fig6_rma0p02_vs_1um_bound` compare all 18 cells against the strict $[0, 1\ \mu\text{m}]^3$ Popoviciu upper bound $R_{g,\max} = \sqrt{3}/2 = 0.8660\ \mu\text{m}$ (`rg_cubic_bound_sweep_eq13`). All **18/18 cells** are **`not reproduced`**, exceeding $0.866\ \mu\text{m}$ by $+0.138$ to $+0.931\ \mu\text{m}$. Combined with Fig. 8 (`cmp_fig6_vs_fig8_unif_1um`, `cmp_fig6_vs_fig8_unif_3um`), the implied initial radius of gyration is $R_{g,i} = 1.532 \pm 0.063\ \mu\text{m}$ across all six $R_{\alpha:a}$ values—$3.06\times$ the uniform $1\ \mu\text{m}^3$ cube value ($0.500\ \mu\text{m}$) and matching the uniform $3\ \mu\text{m}^3$ cube value ($1.500\ \mu\text{m}$) (`q_digitize_fig6_fig8_rg_3x_scaling`).

### Main Fig. 7 (p. 20) — Grid of single trajectory 3D snapshots
- **Fig. 7A–I** *(left out)*: $3 \times 3$ matrix of single trajectory 3D network snapshots at $R_{m:a} \in \{0.005, 0.01, 0.02\}$ and $R_{\alpha:a} \in \{0.01, 0.1, 0.5\}$ after $2000\ \text{s}$. Qualitative snapshots with no numeric axes.

### Main Fig. 8 (p. 21) — Normalized $R_{g,f}/R_{g,i}$ vs time at $R_{m:a} = 0.01$
- **Fig. 8** *(read — 6 readings, 62 points total)*:
  - **`fig8_rg_ratio_ra0p01`**, **`fig8_rg_ratio_ra0p02`**, **`fig8_rg_ratio_ra0p05`**, **`fig8_rg_ratio_ra0p1`**, **`fig8_rg_ratio_ra0p2`**, **`fig8_rg_ratio_ra0p5`** (`fig8_rg_ratio_*_fig`, `run`, `figures/fig8_rg_ratio_*.py` via `_timeseries.py`, 900 dpi, crop `[2380, 800, 7240, 4720]`): All six solid mean trajectories traced using swatch RGBs extracted directly from the figure legend; axes calibrated on 5 x-ticks (`0..2000 s`, residual `0.30 px`) and 8 y-ticks (`0.75..1.10`, residual `0.32 px`).
  - Notes on occlusion and SD bands: `R_alpha:a = 0.01` (red) and `0.02` (green) pass behind the upper-right legend box over $t \in [1550, 1970]\ \text{s}$ before emerging at $t \in [1978, 2000]\ \text{s}$ (`0.9943` and `1.0209`), so those two curves are read at `100, 200, 400, 600, 800, 1000, 1200, 1400, 2000 s` (the other four are read at `100, 200, 400, ..., 2000 s`). All six curves coincide at $(0, 1.00)$ by definition; the six overlapping translucent SD bands cannot be separated uniquely per curve and are omitted from per-curve bounds.
  - Terminal ($t = 2000\ \text{s}$) ratios: `0.9943, 1.0209, 0.8871, 0.8337, 0.8265, 0.8287` for $R_{\alpha:a} = 0.01, 0.02, 0.05, 0.1, 0.2, 0.5$, tested against Fig. 6 in `cmp_fig6_vs_fig8_unif_1um` and `cmp_fig6_vs_fig8_unif_3um`.

### Main Fig. 9 (p. 22) — Heat map of nematic order parameter $S$ after 2000 s
- **Fig. 9** *(read — 3 readings, 18 cells total)*:
  - **`fig9_S_rma0p005`**, **`fig9_S_rma0p01`**, **`fig9_S_rma0p02`** (`fig9_S_rma*_fig`, `run`, `figures/fig9_S_rma*.py` via `_heatmaps.py`, 900 dpi, crop `[2400, 820, 6900, 4320]`): All three rows ($R_{m:a} = 0.005, 0.01, 0.02$) across all six $R_{\alpha:a}$ columns, decoded directly from the image's colorbar calibrated on all 7 printed labels `0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1` (exact uniform spacing `468.0 px` per `0.1`, residual `0.00 px`; bottom label `0.1` aligns with the bottom border while `0.7` sits `72 px` below the top border, giving an exact colorbar span `[0.100, 0.7154]`).
  - Readings:
    - $R_{m:a} = 0.005$: `0.103, 0.103, 0.145, 0.205, 0.239, 0.302`
    - $R_{m:a} = 0.01$: `0.110, 0.106, 0.292, 0.376, 0.473, 0.443`
    - $R_{m:a} = 0.02$: `0.106, 0.239, 0.530, 0.635, 0.713, 0.657`
  - **Cross-figure validation**: At the shared condition ($R_{m:a} = 0.02, R_{\alpha:a} = 0.1, \chi = 1, t = 2000\ \text{s}$), Fig. 9 gives `S = 0.635 [0.619, 0.651]`, agreeing within `0.010` with Fig. 10B (`0.624 ± 0.082` SD).

### Main Fig. 10 (p. 23) — $R_{g,f}/R_{g,i}(t)$ and $S(2000\ \text{s})$ vs turnover factor $\chi$
- **Fig. 10A** *(read — 6 readings, 58 points total)*:
  - **`fig10a_rg_ratio_chi0p125`**, **`fig10a_rg_ratio_chi0p25`**, **`fig10a_rg_ratio_chi0p5`**, **`fig10a_rg_ratio_chi2`**, **`fig10a_rg_ratio_chi4`**, **`fig10a_rg_ratio_chi8`** (`fig10a_rg_ratio_*_fig`, `run`, `figures/fig10a_rg_ratio_*.py` via `_timeseries.py`, 900 dpi, crop `[360, 1020, 3400, 3920]`): All six plotted curves (note that $\chi = 1$ is omitted from Fig. 10A per `q_audit_crossref_and_omissions`); axes calibrated on 5 x-ticks (`0..2000 s`, residual `0.40 px`) and 7 y-ticks (`0.6..1.2`, residual `0.38 px`). Curves $\chi = 0.125, 0.25, 0.5, 2$ are read over `100..2000 s` (11 points each; terminal $R_{g,f}/R_{g,i}$ = `0.6885, 0.6939, 0.7611, 0.8916`); $\chi = 4$ and $8$ enter behind the legend box after `1200 s` and are read over the unoccluded window `100..1200 s` (7 points each). Overlapping translucent SD bands are omitted per curve.
- **Fig. 10B** *(read — 1 reading, 7 points with SD error bars)*:
  - **`fig10b_S_vs_chi`** (`fig10b_S_vs_chi_fig`, `run`, `figures/fig10b_S_vs_chi.py` via `_logx_errbars.py`, 900 dpi, crop `[3400, 1020, 6100, 3920]`, `error_bars: sd`): Log10 x-axis calibrated on 3 major decade ticks (`0.1, 1, 10`, residual `0.17 px`), linear y-axis on 8 ticks (`0.1..0.8`, residual `0.45 px`). Circular marker centers isolated by morphological disk opening and error-bar caps verified symmetric to `<= 0.0002`:
    - $\chi = 0.125$: `0.3757 ± 0.1371`
    - $\chi = 0.25$: `0.3372 ± 0.1065`
    - $\chi = 0.5$: `0.4719 ± 0.1717`
    - $\chi = 1.0$: `0.6244 ± 0.0815`
    - $\chi = 2.0$: `0.4860 ± 0.1397`
    - $\chi = 4.0$: `0.3105 ± 0.1201`
    - $\chi = 8.0$: `0.3049 ± 0.1207`

### Main Fig. 11 (p. 23) — Single trajectory 3D snapshots vs $\chi$
- **Fig. 11A–C** *(left out)*: Three single trajectory 3D snapshots at $\chi = 0.125, 1, 8$ after $2000\ \text{s}$ with black/white plus/minus-end beads. Qualitative snapshots with no numeric axes.

### Main Fig. 12 (p. 24) — Actin filament MSD $\langle\Delta x^2\rangle(t)$ and diffusion exponent $\nu(\chi)$
- **Fig. 12A** *(read — 7 readings, 85 points total)*:
  - **`fig12a_msd_chi0p125`**, **`fig12a_msd_chi0p25`**, **`fig12a_msd_chi0p5`**, **`fig12a_msd_chi1`**, **`fig12a_msd_chi2`**, **`fig12a_msd_chi4`**, **`fig12a_msd_chi8`** (`fig12a_msd_*_fig`, `run`, `figures/fig12a_msd_*.py` via `_timeseries.py`, 900 dpi, crop `[360, 1020, 3080, 3800]`): All seven curves calibrated on 5 x-ticks (`0..2000 s`, residual `0.15 px`) and 6 y-ticks (`0.0..0.5 um^2`, residual `0.20 px`). Caption states SE error bars are smaller than the data points (not visible). Reading windows avoid occlusion: $\chi = 0.125$ (red, occluded by green until ~650 s) read at `700..2000 s` (9 pts); $\chi = 0.25, 0.5$ read at `200..2000 s` (14 pts each); $\chi = 1, 2, 4, 8$ (entering legend box after 1200 s) read at `100..1200 s` (12 pts each).
  - **Model comparisons**: All seven curves are tested against `anomalous_diffusion_msd_eq15` (`prefactor: unknown`, `cmp_eq15_msd_chi*`) using Fig. 12B's exact $\nu(\chi)$ extracted via `extract_nu_at_chi`, and log-log OLS slopes over $t \le 1000\ \text{s}$ (`msd_loglog_fit_eq15`) are compared against Fig. 12B (`cmp_fig12b_nu_vs_fig12a_fit`). Single power-law scaling holds closely up to $1000\ \text{s}$ for $\chi = 0.5$ (`9/9` points in band) and $\chi = 0.125$ (`4/4` points), but departs strongly at $\chi = 4, 8$ where MSD bends over and peaks near $900\ \text{s}$ before the stated $1000\ \text{s}$ cutoff (`q_digitize_fig12_peak_and_powerlaw_fit`).
- **Fig. 12B** *(read — 1 reading, 7 points with regression SE error bars)*:
  - **`fig12b_nu_vs_chi`** (`fig12b_nu_vs_chi_fig`, `fitted`, `figures/fig12b_nu_vs_chi.py` via `_logx_errbars.py`, 900 dpi, crop `[3040, 1020, 5980, 3800]`, `error_bars: sem`): Log10 x-axis calibrated on 3 decade ticks (residual `0.17 px`), linear y-axis on 7 ticks (`0.9..1.5`, residual `0.21 px`). Error-bar cap asymmetry `<= 0.0001` across all visible bars; at $\chi = 0.125$, the error bar is hidden inside the circular marker and reported as at most the marker radius (`<= 0.0097`):
    - $\chi = 0.125$: `1.1235` (bar hidden inside marker, `<= 0.0097`)
    - $\chi = 0.25$: `1.2778 ± 0.0181`
    - $\chi = 0.5$: **`1.4107 ± 0.0196`** *(true peak, above $\chi = 1$, vs text p. 24 "centered around $\chi = 1$")*
    - $\chi = 1.0$: `1.3805 ± 0.0266`
    - $\chi = 2.0$: `1.3366 ± 0.0462`
    - $\chi = 4.0$: `1.1686 ± 0.0504`
    - $\chi = 8.0$: `1.0503 ± 0.0624`

### Main Fig. 13 (p. 26) — Large $3 \times 3 \times 3\ \mu\text{m}^3$ network 3D snapshot
- **Fig. 13** *(left out)*: Single trajectory 3D snapshot of a $3 \times 3 \times 3\ \mu\text{m}^3$ actomyosin system at $R_{\alpha:a} = 0.1, R_{m:a} = 0.02$ after $500\ \text{s}$, color-coded by XY-angle ($-90^\circ$ to $90^\circ$). Qualitative spatial rendering without numeric plot axes or aggregated angle distribution.

---

## Supporting Information Figures (PDF pp. 29, 39, 50–54)

### External PLOS SI Figures & Videos (captions on p. 29; files not bundled in PDF)
- **S1 Fig** *(left out as external PNG not bundled in the 62-page PDF)*: Three 3D snapshots illustrating $R_g = 1.77, 1.51, 1.22\ \mu\text{m}$ (printed caption values already recorded as `s1fig_rg_a`, `s1fig_rg_b`, `s1fig_rg_c` and compared against the $0.866\ \mu\text{m}$ Popoviciu bound).
- **S2 Fig** *(left out as external PNG not bundled in the 62-page PDF)*: Binding-site exclusivity heat map and $R_{g,f}/R_{g,i}(t)$ time series; contradiction between Main Text p. 25 and S2 Fig caption p. 29 is logged in `q_audit_s2fig_exclusivity_contradiction`.
- **S1–S4 Videos** *(left out as external MOV files not bundled in the PDF)*.

### S1 Text Fig. 1 (p. 39 / S1-4) — Software architecture workflow
- **S1 Text Fig. 1** *(left out)*: Block workflow diagram describing the C++ MEDYAN software architecture (User input, System, Structure, Controller). No quantitative data.

### S4 Text Fig. 1 (p. 50 / S4-2) — Mean filament length vs simulation time (NRM vs Gillespie DM)
- **S4 Fig. 1** *(read — 2 readings, 22 points total)*:
  - **`s4fig1_length_nrm`**, **`s4fig1_length_gillespie`** (`s4fig1_length_*_fig`, `run`, `figures/s4fig1_length_*.py` via `_s4_excel.py`, 900 dpi, crop `[1320, 2760, 6280, 5960]`): Calibrated from 10 horizontal Excel gridlines (`0.0..0.9 um`, residual `1.47 px`) and 11 x-axis labels (`0..10 s`). No error bars are plotted. Both algorithms rise from `0.0275 um` at $t = 0\ \text{s}$ to **`0.7899 um`** at $t = 10\ \text{s}$ (NRM slightly below Gillespie at $t = 1\text{–}4\ \text{s}$: `0.4688, 0.6011, 0.6663, 0.7102` vs `0.4786, 0.6126, 0.6777, 0.7177 um`, and coincident within marker radius for $t = 5\text{–}10\ \text{s}$).
  - **Model comparisons**: `cmp_s4fig1_nrm_final_vs_printed` and `cmp_s4fig1_gillespie_final_vs_printed` confirm the $t = 10\ \text{s}$ terminal length (`0.790 um`) reproduces Main Text p. 17's `0.8 um`, confirming that the ~19% excess over strict monomer mass balance (`0.634–0.644 um`) is already present in the pure-actin 10 s benchmark (`q_digitize_s4fig1_s4fig2_equilibration_and_mass_balance`).

### S4 Text Fig. 2 (p. 51 / S4-3) — Diffusing G-actin concentration vs simulation time (NRM vs Gillespie DM)
- **S4 Fig. 2** *(read — 2 readings, 22 points total with unstated error bars)*:
  - **`s4fig2_conc_nrm`**, **`s4fig2_conc_gillespie`** (`s4fig2_conc_*_fig`, `run`, `figures/s4fig2_conc_*.py` via `_s4_excel.py`, 900 dpi, crop `[1320, 940, 6280, 4580]`, `error_bars: unstated`): Calibrated from 6 horizontal Excel gridlines (`0..25 uM`, residual `0.24 px`) and 11 x-axis labels (`0..10 s`). Both series decay from `19.32 uM` at $t = 0\ \text{s}$ (`8.06 ± 0.32` NRM / `8.28 ± 0.44` Gillespie at 1 s; `4.73 ± 0.37` / `5.00 ± 0.46` at 2 s; `3.10 ± 0.34` / `3.32 ± 0.33` at 3 s) to **`0.492 ± 0.14 uM`** (NRM) and **`0.492 ± 0.12 uM`** (Gillespie) at $t = 10\ \text{s}$ (with error-bar ends occluded by the overlapping series mirrored symmetrically as documented in each point's note; error bars at $t = 0$ and small late bars hidden inside the marker reported as at most marker radius).
  - **Model comparisons**: `cmp_s4fig2_nrm_final_vs_analytical_Ac`, `cmp_s4fig2_gillespie_final_vs_analytical_Ac` (`not reproduced`: `0.492 uM` vs analytical `0.17 uM`), and `cmp_s4fig2_nrm_final_vs_sim_claim` (`0.492 uM` vs `"about 0.4 uM"`, still declining from `0.693 uM` at 8 s and `0.559 uM` at 9 s), quantitatively settling `q_audit_ss_actin_conc_s2_vs_s4` (`q_digitize_s4fig1_s4fig2_equilibration_and_mass_balance`).

### S4 Text Fig. 3 (p. 52 / S4-4) — Reaction-diffusion wall-clock computation time (NRM vs Gillespie DM)
- **S4 Fig. 3** *(read — 2 readings, 22 points total; Gillespie with unstated error bars)*:
  - **`s4fig3_time_nrm`**, **`s4fig3_time_gillespie`** (`s4fig3_time_*_fig`, `run`, `figures/s4fig3_time_*.py` via `_s4_excel.py`, 900 dpi, crop `[1320, 940, 6280, 4580]`): Calibrated from 8 horizontal Excel gridlines (`0..7000 s`, residual `1.00 px`) and 11 x-axis labels (`0..10 s`).
  - Readings:
    - Gillespie DM (red, `error_bars: unstated`): `0.0, 998.8, 2531.6, 3662.0 ± 154, 4377.5 ± 174, 4887.6 ± 229, 5255.4 ± 257, 5521.1 ± 255, 5725.1 ± 225, 5878.2 ± 210, 6000.7 ± 218 s`.
    - NRM (blue, no visible error bars): visibly elevated above the $y = 0$ gridline across $t = 1..10\ \text{s}$: `0.0, 29.0, 49.2, 69.6, 69.6, 79.9, 90.2, 90.0, 100.3, 100.2, 100.3 s`.
  - **Key quantitative finding**: `cmp_s4fig3_speedup_t10_vs_claim` and `cmp_s4fig3_speedup_max_vs_claim` yield a Gillespie-to-NRM wall-clock speedup of **`59.8-fold`** at $t = 10\ \text{s}$ (`6000.7 / 100.3 s`) and a peak of **`62.9-fold`** at $t = 4\ \text{s}$ (`4377.5 / 69.6 s`) — well below the `"over 100-fold"` performance increase claimed in both S4 Text Section A (p. 50) and the S4 Fig. 3 caption (p. 52) (`q_digitize_s4fig3_speedup_overstated`).

### S4 Text Fig. 4 (p. 53 / S4-5) — Polymer midpoint displacement error $e_d$ vs segment length
- **S4 Fig. 4** *(read — 3 readings, 12 points total)*:
  - **`s4fig4_err_actin_5pN`**, **`s4fig4_err_actin_1pN`**, **`s4fig4_err_mt_1pN`** (`s4fig4_err_*_fig`, `run`, `figures/s4fig4_err_*.py` via `_s4_excel.py`, 900 dpi, crop `[1320, 940, 6280, 4580]`): Calibrated from 8 horizontal Excel gridlines (`0.00..0.07`, residual `1.07 px`) across the 4 plotted coarse-graining levels (`1` [Bead-spring], `10`, `20`, `40` monomers). No error bars are plotted.
  - Readings:
    - `Actin filament, 5 pN tolerance` (red): `0.0000, 0.04872 (4.87%), 0.06079 (6.08%), 0.06638 (6.64%)`
    - `Actin filament, 1 pN tolerance` (blue): `0.0000, 0.00922 (0.92%), 0.03577 (3.58%), 0.04323 (4.32%)`
    - `Microtubule, 1 pN tolerance` (green): `0.0000, 0.00520 (0.52%), 0.02608 (2.61%), 0.03951 (3.95%)`
  - **Key quantitative finding**: `cmp_s4fig4_peak_actin_5pN_vs_claim`, `cmp_s4fig4_peak_actin_1pN_vs_claim`, and `cmp_s4fig4_peak_mt_1pN_vs_claim` quantitatively settle `q_s4_fig4_displacement_error_overstated`: `Actin 5 pN` peaks at **`6.64%`** (`+1.64%` above `5%`, with `6.08%` already above `5%` at 20 monomers), directly refuting the caption claim `"Displacement errors stay under 5% of the total polymer length for all levels of coarse-graining used"`, whereas both `1 pN` series remain below `5%` (`q_digitize_s4fig4_error_exceeds_5pct`).

### S4 Text Fig. 5 (p. 54 / S4-6) — Single-filament equilibration wall-clock time vs segment length
- **S4 Fig. 5** *(read — 3 readings, 12 points total)*:
  - **`s4fig5_eqtime_actin_5pN`**, **`s4fig5_eqtime_actin_1pN`**, **`s4fig5_eqtime_mt_1pN`** (`s4fig5_eqtime_*_fig`, `run`, `figures/s4fig5_eqtime_*.py` via `_s4_excel.py`, 900 dpi, crop `[1320, 940, 6280, 4580]`): Calibrated from 7 horizontal Excel gridlines (`0..6 s`, residual `1.64 px`) across `1, 10, 20, 40` monomers. No error bars are plotted.
  - Readings:
    - `Actin filament, 5 pN tolerance` (red): `0.2217, 0.1497, 0.1276, 0.0428 s`
    - `Actin filament, 1 pN tolerance` (blue): `3.1373, 1.7333, 0.3502, 0.1414 s`
    - `Microtubule, 1 pN tolerance` (green): `5.3437, 3.8543, 1.9782, 0.2095 s`
