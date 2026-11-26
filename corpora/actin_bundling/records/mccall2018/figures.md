# Figure Digitization Log — McCall et al. (2018)

Source PDF: `corpus/actin_bundling/2018_McCall.pdf` (33 pages: Main Text *Biophysical Journal* 114, 1636–1645, PDF pp. 1–10; Supporting Material cover and SI-1 through SI-22, PDF pp. 11–33).

- **Summary:** Digitized **42 calibrated figure readings (652 data points)** spanning all **17 quantitative plot subpanels** across Main Text Figs. 1, 3, 4 and Supporting Figs. S1–S6. Main Text panels are crisp vector PDF graphics rendered at **900–1200 dpi**; Supporting Material panels are embedded ~150–161 ppi JPEGs rendered at **450–600 dpi**. Every reading uses self-contained image-based frame/tick calibration and curve/marker/bar detection in `figures/_calib.py` + `figures/<id>.py` (no PDF re-rendering inside `read(image, at)`), passed `boxgraph figure check` with 0 errors, and was visually verified point-by-point on compact overlay contact sheets in `/tmp/figwork_mccall2018/sheets/`. Saved frozen readings in `readings.json`.

---

## Main Text Figures (Figs. 1–4)

### Figure 1 (PDF p. 2 / printed p. 1637)
- **Fig. 1A, 1B, 1C:** Left out — spinning-disk confocal fluorescence (TMR-actin, 647-BSA) and DIC micrographs (coverslip surface, droplet midplane, and x-z cross section; qualitative imaging with no quantitative plot axes).
- **Fig. 1D (READ — 3 readings, 63 pts; 1200 dpi):** `fig1d_actin` (21 pts, green curve), `fig1d_bsa` (21 pts, orange curve), `fig1d_dic` (21 pts, black curve) across `-5.5` to `+5.5 µm` (`I_norm` in `[0, 1]`). Calibrated from frame spines and inward bottom/left ticks (`cal_res < 0.4 px`). Captures twin peripheral actin peaks (`~0.92` at `±2.75 µm`, interior `~0.08`) vs uniform BSA interior plateau (`~0.85–0.96`).
- **Fig. 1E (READ — 2 readings, 4 bars + SEM error bars; 1200 dpi):** `fig1e_actin_ram` (green bars: `x=1` TMR-actin + BSA `n=171` = `33.23 ± 0.28`, `x=2` TMR-actin alone `n=117` = `36.45 ± 0.57`) and `fig1e_bsa_ram` (orange bars: `x=3` 647-BSA + actin `n=214` = `8.88 ± 0.11`, `x=4` 647-BSA alone `n=287` = `8.76 ± 0.13`). Calibrated from right major ticks `0..40` (`cal_res = 0.40 px`); bar tops and SEM caps verified on zoomed overlays. **Key finding (`q_fig1e_actin_ram_alone_vs_bsa_and_30_rounding`):** BSA bars agree within `0.7 SEM` (`~10%` above printed `"≈ 8"`), whereas actin alone (`0.5 µM`) is `9.7%` (`5.0 combined SEMs`) higher than actin + BSA (`0.25 µM` actin + `0.25 µM` BSA), challenging the p. 3 statement that actin RAM is `"the same whether one or both proteins are present"`. Actin/BSA ratios are `3.74×` and `4.16×` (mean `3.95× ≈ 4`).

### Figure 2 (PDF p. 3 / printed p. 1638)
- **Fig. 2A, 2B, 2C, 2D:** Left out — confocal fluorescence (TMR-actin, 647-phalloidin, merge) and DIC micrographs (qualitative imaging).
- **Fig. 2E:** Left out — 2D pixel-cloud count-density heatmap (`I_647-Phalloidin` vs `I_TMR-actin`, color scale 0–60 counts; thousands of overlapping pixel bins with no discrete points or fitted curve; summary Pearson correlation `r = 0.864` is already recorded as scalar datum `r_actin_phalloidin`).

### Figure 3 (PDF p. 4 / printed p. 1639)
- **Fig. 3A:** Left out — conceptual cartoon schematic of spontaneous pyrene-actin polymerization phases and `t_1/2` definition (no numeric axes).
- **Fig. 3B (READ — 3 readings, 48 pts; 900 dpi):** `fig3b_actin_alone` (16 pts, black: `3.35 -> 13.47 a.u.`), `fig3b_plus_plk` (16 pts, cyan: `6.65 -> 35.21 a.u.`), `fig3b_plus_coacervate` (16 pts, red: `7.94 -> 34.56 a.u.`) over `1–29.5 min`. Calibrated from top/left major ticks (`cal_res = (0.21, 0.00) px`). **Key finding (`q_fig3b_vs_fig3c_halftime_discrepancy`):** applying SI Eq. (5) (`I_1/2 = (I_max + I_0)/2`, `I_0 = 3.35`) gives `1/t_1/2 = 0.312 min^-1` (`t_1/2 = 3.21 min`) for `+3 µM pLK+pRE` and `0.230 min^-1` (`4.35 min`) for `+3 µM pLK`, `45%–53%` faster than the `3.0 µM` points plotted in Fig. 3C (`0.215` and `0.150 min^-1`).
- **Fig. 3C (READ — 4 readings, 30 pts; 900 dpi):** Log-log calibration from decade major ticks and minor ticks (`cal_res = (0.35, 0.36) px`). `fig3c_coacervate` (13 red circles, `0.3–3000 µM`: `0.0401 -> 2.854 min^-1`), `fig3c_plk_alone` (11 cyan squares, `0.3–30 µM`: `0.0447 -> 0.601 min^-1` at `15 µM`, `0.546 min^-1` at `30 µM`, including red-occluded squares cleanly centered via edge/4-corner geometry), `fig3c_control_coacervate` (3 pts on red dashed baseline `0.0281 min^-1`), `fig3c_control_plk` (3 pts on cyan dashed baseline `0.0249 min^-1`). **Key finding (`q_fig3c_saturation_above_30uM_overstated`):** rate increases from `1.20 min^-1` (`42.9×` control) at `30 µM` to `1.44 min^-1` (`51.1×`) at `300 µM` and `2.85 min^-1` (`101.7×`) at `3000 µM` — a further `2.37×` rise rather than flat saturation above `30 µM`.
- **Fig. 3D (READ — 3 readings, 48 pts; 900 dpi):** `fig3d_actin_alone` (16 pts, black: flat `11.80–12.84 a.u.`), `fig3d_plus_plk` (16 pts, cyan: rapid decay from `33.41` at `0.5 min` to `23–24 a.u.`), `fig3d_plus_coacervate` (16 pts, red: mild rise `14.5 -> 19.1 a.u.`). Calibrated from top/right major ticks (`cal_res = (0.64, 0.41) px`).

### Figure 4 (PDF p. 5 / printed p. 1640)
- **Fig. 4A, 4B, 4C:** Left out — confocal fluorescence micrographs (surface, midplane, magnified insets) across `[TMR-actin] = 0.01–1.5 µM` (qualitative imaging).
- **Fig. 4D (READ — 6 readings, 114 pts; 1200 dpi):** Semi-log intensity profiles (`-3.8` to `+3.8 µm`, `10^0–4×10^4 counts`) separated cleanly by distinct RGB hues: `fig4d_1p5uM` (purple), `fig4d_0p5uM` (blue), `fig4d_0p25uM` (green), `fig4d_0p10uM` (yellow), `fig4d_0p05uM` (orange), `fig4d_0p01uM` (red), 19 pts each. Calibrated from top linear and right log10 ticks (`cal_res < 0.5 px`).
- **Fig. 4E (READ — 3 readings, 18 pts; 1200 dpi):** Mean droplet intensity (`<Intensity>_droplet`, counts) vs `[TMR-actin]` (`0.01, 0.05, 0.1, 0.25, 0.5, 1.5 µM`): `fig4e_periphery` (6 cyan triangles: `56.9 -> 20,875`), `fig4e_interior` (6 black circles: `24.6 -> 4,080`), `fig4e_exterior` (6 red squares: `11.5 -> 466`). Calibrated from top/right log10 major ticks (`cal_res < 0.4 px`).
- **Fig. 4F (READ — 3 readings, 19 pts; 1200 dpi):** `fig4f_ram_periphery` (6 cyan triangles: `5.04, 16.19, 42.10, 46.23, 46.23, 43.76`), `fig4f_ram_interior` (6 black circles: `2.16, 3.86, 5.32, 5.69, 7.44, 8.62`), and `fig4f_filament_boundary` (7 pts along red dashed line, `0.01–0.09 µM`: `10.005 -> 1.111`). **Model comparisons:** `cmp_fig4f_filament_boundary` reproduces `filament_boundary_ram` (`RAM_crit = 0.1 / c_global`) within `0.02%–0.05%` across all 7 points (`reproduced`); `cmp_fig4ef_aux_ram_periphery` and `cmp_fig4ef_aux_ram_interior` confirm internal consistency of `auxiliary_ram_ratio` (`I_region / I_exterior` from Fig. 4E vs Fig. 4F) to `< 0.1%–2.3%` across all 12 points (`close`).

---

## Supporting Material Figures (Figs. S1–S6 & Movies S1–S2)

### Figure S1 (PDF p. 25 / printed p. SI-14)
- **Fig. S1A:** Left out — DIC coalescence time-lapse micrographs (`-10` to `60 ms`; merged diameter `d = 4 µm` already recorded as scalar datum `d_coalesced_figS1`).
- **Fig. S1B (READ — 2 readings, 17 pts + dy error bars; 600 dpi):** `figS1b_aspect_ratio_data` (7 black squares at `0.00–0.06 s` with vertical `dy` uncertainty bars: `1.676 ± 0.083, 1.153 ± 0.065, 1.003 ± 0.060, 1.018 ± 0.061, 0.972 ± 0.060, 0.992 ± 0.060, 1.014 ± 0.061`; note horizontal `dx` acquisition-time bars are also drawn) and `figS1b_aspect_ratio_fit` (10 pts along red exponential curve, `1.663 -> 1.001`). Calibrated from all 7 top x-ticks and 6 left y-ticks (`cal_res = (0.41, 0.51) px`). **Model comparisons:** `cmp_figS1b_exponential_fit` reproduces `1 + (AR_0 - 1)*exp(-t / 6.33 ms)` within `±0.0002–0.0072` across all 10 points (`reproduced`), and `cmp_figS1b_data_vs_fit` confirms all 7 measured points agree with the fit well inside their `dy` error bars (`reproduced`).

### Figure S2 (PDF p. 26 / printed p. SI-15)
- **Fig. S2A, S2B:** Left out — confocal fluorescence and DIC micrographs for actin alone and BSA alone (qualitative imaging).
- **Fig. S2C (READ — 4 readings across 2 stacked subpanels, 60 pts; 600 dpi):** Top subpanel (`figS2c_actin`, 15 pts green; `figS2c_actin_dic`, 15 pts black) and bottom subpanel (`figS2c_bsa`, 15 pts orange; `figS2c_bsa_dic`, 15 pts black) across `-3.6` to `+3.6 µm`. Calibrated from the 3 horizontal spines, 2 vertical spines, and inward top/left/right ticks (`cal_res < 0.7 px`).

### Figure S3 (PDF pp. 27–28 / printed pp. SI-16 to SI-17)
- **Fig. S3A, S3C, S3E, S3F:** Left out — flow-in confocal fluorescence and DIC micrographs and maximum-intensity z-projection (qualitative imaging).
- **Fig. S3B (READ — 1 reading, 20 pts; 600 dpi):** `figS3b_bsa_flowin_ram` (all 20 circular markers: 5 pre-flow points `-8..0 s` at `RAM = 0.90–0.91`, 15 post-flow points `35..63 s` rising smoothly from `6.52` to `9.04`). Calibrated from frame spines and bottom/right ticks (`cal_res < 0.6 px`).
- **Fig. S3D (READ — 2 readings, 24 pts; 600 dpi):** `figS3d_actin_flowin_profile` (21 pts along black intensity profile across `-14` to `+14 µm`: exterior `~106–115`, twin peripheral peaks `~497` and `~593` counts near `±6.2 µm`, interior `~148–160` counts) and `figS3d_exterior_baseline` (3 pts on red dashed exterior guide line at `110.4` counts). Calibrated from bottom/right major ticks (`cal_res < 0.6 px`).

### Figure S4 (PDF p. 29 / printed p. SI-18)
- **Fig. S4A, S4B:** Left out — confocal fluorescence and DIC micrographs for post-phase-separation TMR-actin addition (qualitative imaging).
- **Fig. S4C (READ — 2 readings, 30 pts; 600 dpi):** `figS4c_actin` (15 pts, green) and `figS4c_dic` (15 pts, black) across `-9` to `+9 µm`. Calibrated from frame spines and inward ticks (`cal_res < 0.7 px`).
- **Fig. S4D (READ — 1 reading, 95 resolved scatter points; 600 dpi):** `figS4d_ram_vs_radius` (`RAM` vs `Droplet Radius` across `1–5 µm`). Distance-transform peak detection resolves 95 distinct circular dot centers out of `N = 117` reported droplets (~22 partially overlapping pairs in the dense `R = 1.5–3 µm` core; 0 false positives verified visually). **Settles `q_figs4d_ram_magnitude_and_sem_vs_sd` (`q_figs4d_digitized_confirms_sd_not_sem`):** digitized sample mean = `70.26` (printed `70.3`), sample **SD = `16.59`** (printed `15.9`), sample **SEM = `1.70`** (`1.47` for `N = 117`), and Pearson `r = -0.448` (printed `-0.427`), confirming `± 15.9` is SD rather than SEM.

### Figure S5 (PDF p. 30 / printed p. SI-19)
- **Fig. S5A, S5B, S5C, S5D:** Left out — FRAP vs control DIC/fluorescence micrographs and focal-plane schematic (qualitative imaging).
- **Fig. S5E (READ — 2 readings, 40 pts; 450 dpi):** `figS5e_frap_actin` (20 pts, green squares, `-8..238 s`: `1.000 -> ~0.00` with `150–238 s` mean `0.031`, max `0.084`) and `figS5e_frap_bsa` (20 pts, orange circles, `-8..238 s`: `1.000 -> 0.011`, recovering to `0.287` at `30 s`, `0.420` at `90 s`, jumping to `0.533` at `105 s`, and plateauing at `0.604` mean over `150–238 s`). Calibrated from all 9 top major x-ticks and 6 right major y-ticks (`cal_res < 0.5 px`).

### Figure S6 (PDF p. 31 / printed p. SI-20)
- **Fig. S6A, S6B, S6D:** Left out — fluorescence time-lapse micrographs of actin filament winding (± methylcellulose) and circular-arc kymograph (qualitative imaging).
- **Fig. S6C (READ — 1 reading, 22 pts; 600 dpi):** `figS6c_phalloidin_winding` (22 pts at `100..2200` along black step-like trace, `3.86 -> 10.70 × 10^2` counts, capturing baseline `~3.87` and successive winding plateaus `~4.7, ~6.1, ~7.6, ~9.2, ~10.7 × 10^2`). X-axis unit recorded verbatim as printed (`min`), noting the evident misprint for seconds (`s`) per `q_figs6c_time_unit_misprint`. Calibrated from bottom/right major ticks (`cal_res < 0.6 px`).

### Supplemental Movies S1 & S2 (PDF pp. 32–33 / printed pp. SI-21 to SI-22)
- **Movie S1, Movie S2 stills:** Left out — single video-frame preview thumbnails (qualitative).
