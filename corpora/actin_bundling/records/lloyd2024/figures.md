# Figure Inventory and Digitization Log — `lloyd2024`

Complete panel-by-panel inventory of `corpus/actin_bundling/2024_Lloyd.pdf` (37 PDF pages: pp. 1–15 Nature Communications 15, 4073 (2024) main article & Methods; pp. 16–37 Supplementary Information containing Supplementary Results, Supplementary Methods, Supplementary Figs. 1–12, and Supplementary References).

All **45 calibrated figure readings** in [`figures.yaml`](figures.yaml) were implemented as 900–1800 dpi `mode: code` Python readers under [`figures/`](figures/), verified with `boxgraph figure check`, visually inspected on every generated magenta overlay PNG in `figwork/`, persisted in [`readings.json`](readings.json), linked to `origin: measured` entries in [`data.yaml`](data.yaml), and exercised by the record's curve-port and scaling models in [`runs.yaml`](runs.yaml) and [`comparisons.yaml`](comparisons.yaml).

---

## Page 1 (printed p. 1)
- No figures (Title, Abstract, Introduction text only).

## Page 2 (printed p. 2) — Main Figure 1: The neotropical butterfly *Heliconius sara*
- **Fig. 1A–C**: Left out — dorsal photographs of *H. sara* and magnified iridescent blue / non-iridescent black-and-white wing regions (no plot axes).
- **Fig. 1D–F**: Left out — scanning electron micrographs (SEMs) of overlapping cover/ground scales, longitudinal ridges, crossribs, windows, and microribs (scale bars only).
- **Fig. 1G**: Left out — transmission electron micrograph (TEM) cross-section through scale ridges showing layered lamellae (no plot axes).

## Page 3 (printed p. 3)
- No figures (Results text only).

## Page 4 (caption on Page 5) — Main Figure 2: Adult ridge morphology and 50% pupal actin patterning
- **Fig. 2A–C**: Left out — false-colored and representative SEM micrographs illustrating scale length, width, area, ridge spacing, and ridge width definitions.
- **Fig. 2D (Adult *H. sara* ridge number vs scale width)**: Digitized both OLS linear regression lines and shaded 95% CI ribbons at 1200 dpi (`rx = 1.60 px`, `ry = 1.20 px`, `error_bars: ci95`; individual semi-transparent scatter points overlap densely within cover/ground clusters, and marginal density curves have no numerical y-axis scale):
  - `fig2d_reg_blue` (`d_fig2d_reg_blue`, 12 points `W = 26..48 µm`, line residual `0.77 px`): rises linearly from `34.50` (`26 µm`) to `58.78` (`48 µm`); at adult blue cover mean `W = 29.5 µm`, regression gives `38.36` ridges (matching Fig. 2L blue adult dashed line `38.61`). Sweeping `ridge_number_from_width_and_spacing` (`N = W / 0.804 µm`, `run_fig2d_ridge_num_blue_sweep`) gives `cmp_fig2d_ridge_num_blue`: **`reproduced`**.
  - `fig2d_reg_black` (`d_fig2d_reg_black`, 13 points `W = 26..50 µm`, line residual `0.88 px`): rises linearly from `26.83` (`26 µm`) to `52.85` (`50 µm`); at adult black cover mean `W = 31.3 µm`, regression gives `32.57` ridges (matching `run_ridge_number_black_cover = 32.54` and Fig. 2L black adult dashed line `32.69`). Sweeping `ridge_number_from_width_and_spacing` (`N = W / 0.962 µm`, `run_fig2d_ridge_num_black_sweep`) gives `cmp_fig2d_ridge_num_black`: **`reproduced`**.
- **Fig. 2E (Ridge number vs ridge spacing)**: Left out — 80 heavily overlapping semi-transparent individual-mean circles and triangles with no fitted regression curve (several concentric markers cannot be separated unambiguously), and marginal density curves have no numerical density axis. The underlying group distributions are captured quantitatively in SI Fig. 5A and Fig. 2D.
- **Fig. 2F (Ridge spacing vs ridge width)**: Digitized both OLS linear regression lines and shaded 95% CI ribbons at 1200 dpi (`rx = 0.00 px`, `ry = 2.00 px`, `error_bars: ci95`):
  - `fig2f_reg_blue` (`d_fig2f_reg_blue`, 6 points `0.27..0.37 µm`, line residual `0.70 px`): rises steeply from `0.684 µm` to `0.956 µm`.
  - `fig2f_reg_black` (`d_fig2f_reg_black`, 7 points `0.275..0.330 µm`, line residual `0.67 px`): rises gently from `0.933 µm` to `0.992 µm`.
- **Fig. 2G–I, Ii**: Left out — confocal fluorescence micrographs of phalloidin-stained 50% pupal *H. sara* forewings and actin bundles.
- **Fig. 2J (50% developing F-actin bundle width notched boxplots)**: Digitized at 1200 dpi (`ry = 0.80 px`, `error_bars: range`) as `fig2j_actin_width_box` (`d_fig2j_actin_width_box`): `1=black` median `0.4526 µm` (range `[0.4342..0.4924 µm]`), `2=blue` median `0.4429 µm` (range `[0.3872..0.4675 µm]`).
- **Fig. 2K (50% developing F-actin bundle spacing notched boxplots)**: Digitized at 1200 dpi (`ry = 0.80 px`, `error_bars: range`) as `fig2k_actin_spacing_box` (`d_fig2k_actin_spacing_box`): `1=black` median `1.121 µm` (range `[1.004..1.575 µm]`), `2=blue` median `1.051 µm` (range `[0.8204..1.319 µm]`).
- **Fig. 2L (50% developing F-actin bundle number notched boxplots & adult cover dashed lines)**: Digitized at 1200 dpi (`ry = 0.00 px`):
  - `fig2l_actin_number_box` (`d_fig2l_actin_number_box`, `error_bars: range`): `1=black` median `31.97` (range `[24.81..42.14]`, integer bundle counts `25..42`), `2=blue` median `39.00` (range `[30.89..54.01]`, integer bundle counts `31..54`).
  - `fig2l_adult_ridge_number_lines` (`d_fig2l_adult_ridge_number_lines`): horizontal dashed adult cover ridge-number reference lines at **`32.69`** (`1=black`, matching `run_ridge_number_black_cover = 32.54`) and **`38.61`** (`2=blue`, matching Fig. 2D blue regression `38.36` and `29.5 / 0.7483 = 39.4`, settling part 2 of `q_actin_vs_adult_ridge_spacing_offset` in `q_digitize_actin_adult_cover_settlement`).

## Page 5 (printed p. 5)
- No new figures (Fig. 2 caption and Results text).

## Page 6 (printed p. 6) — Main Figure 3: Confocal series of normal wing scale development in *H. sara*
- **Fig. 3A–L, Li**: Left out — confocal fluorescence micrographs (DAPI, phalloidin, anti-α-tubulin, WGA) at 25%–56%+ development (no quantitative plot axes).

## Page 7 (printed p. 7) — Main Figure 4: TauSTED microscopy of actin rearrangement in *H. sara*
- **Fig. 4A–L, Ai, Di, Gi, Ji**: Left out — depth-coded and merged TauSTED super-resolution fluorescence micrographs at 44%, 50%, 63%, and 81% development (no quantitative plot axes).

## Page 8 (printed p. 8) — Main Figure 5: TauSTED imaging of actin associated with scale ultrastructures
- **Fig. 5A–C, Ai, Bi, Ci**: Left out — merged TauSTED fluorescence micrographs and qualitative cartoon insets at 63%, 69%, and 75% development (no quantitative plot axes).

## Page 9 (printed p. 9) — Main Figure 6: Actin patterning in *Parides arcas* and *Morpho helenor*
- **Fig. 6A–L, Ai, Di, Gi, Ji**: Left out — depth-coded and merged TauSTED fluorescence micrographs at ~50%–65% development (no quantitative plot axes).

## Page 10 (printed p. 10) — Main Figure 7: TauSTED actin–chitin colocalization in *Parides arcas* (~50%)
- **Fig. 7A–D**: Left out — merged XY, XZ, and YZ TauSTED fluorescence micrographs (no quantitative plot axes).
- **Fig. 7E (XZ intensity profile 1 of actin & chitin, `0.0–6.8 µm`)**: Digitized both curves on a `0.1 µm` grid (69 points each) at 1200 dpi (white-on-black plot calibrated from x-digit and y-digit centroids aligned to the green `y = 0` baseline; `rx = 1.40 px`, `ry = 1.00 px`; baseline regions where green overlays magenta are recovered at `y = 0`):
  - `fig7e_actin` (`d_fig7e_actin`, green phalloidin): peaks at `1.3 µm` (`47.9`), `6.0 µm` (`33.6`), and `5.4 µm` (`24.0`).
  - `fig7e_chitin` (`d_fig7e_chitin`, magenta CBD): peaks at `5.0 µm` (`64.3`), `0.2 µm` (`51.0`), and `2.1 µm` (`47.7`), straddling the actin peaks.
- **Fig. 7F (XZ intensity profile 2 of actin & chitin, `0.0–8.2 µm`)**: Digitized both curves on a `0.1 µm` grid (83 points each) at 1200 dpi (`rx = 0.80 px`, `ry = 1.20 px`):
  - `fig7f_actin` (`d_fig7f_actin`, green phalloidin): peaks at `6.9 µm` (`53.1`), `3.5 µm` (`34.5`), `7.3 µm` (`30.8`), and `0.9 µm` (`20.3`), sitting between each pair of chitin ridge walls.
  - `fig7f_chitin` (`d_fig7f_chitin`, magenta CBD): paired ridge-wall peaks at `(0.6, 1.4)`, `(3.2, 4.6)`, and `(6.6, 7.8) µm` (heights `28–40`).

## Page 11 (printed p. 11) — Main Figure 8: Cytochalasin D perturbation of actin at 50% development
- **Fig. 8A–D**: Left out — dorsal forewing photographs comparing cyto-D-treated and DMSO-control wings (no plot axes).
- **Fig. 8E (Mean ± SE reflectance spectra, `350–700 nm`)**: Digitized both mean curves and shaded `±SE` ribbons every `25 nm` (15 points each) at 1200 dpi (`rx = 1.28 px`, `ry = 1.60 px`, `error_bars: sem`; symmetric mirroring across purple ribbon overlap at `350–425 nm` and `650–700 nm`):
  - `fig8e_reflectance_control` (`d_fig8e_reflectance_control`, red, `n = 15`): rises from `4.20%` (`350 nm`) to a broad blue-green peak of **`11.02%`** (`525 nm`) / `11.01%` (`550 nm`), falling to `6.37%` (`700 nm`).
  - `fig8e_reflectance_cytod` (`d_fig8e_reflectance_cytod`, blue, `n = 21`): flattened spectrum from `5.15%` (`350 nm`) peaking weakly at **`7.31%`** (`575 nm`) and remaining `7.04%` at `700 nm`.
- **Fig. 8F–M**: Left out — SEM micrographs of cyto-D-treated and control scales at three magnifications.
- **Fig. 8N (Ridge number dot plot & red mean ± SD)**: Digitized at 1200 dpi (`ry = 1.40 px`):
  - `fig8n_ridge_number_mean_sd` (`d_fig8n_ridge_number_mean_sd`, `error_bars: sd`): `1=control` **`45.17 ± 6.54`**, `2=cyto-D` **`46.19 ± 8.74`**.
  - `fig8n_ridge_number_control_scales` (`d_fig8n_ridge_number_control_scales`): all **20/20** individual control black dots (`28.99..58.02`, sample mean `45.28`, SD `6.48`).
  - `fig8n_ridge_number_cytod_scales` (`d_fig8n_ridge_number_cytod_scales`): all **20/20** individual cyto-D black dots (`32.49..67.50`, sample mean `46.20`, SD `8.72`). Exact 20 vs 20 dot count quantitatively settles `q_fig8_NO_swapped_refs_and_n` (caption misprints `n = 12`).
- **Fig. 8O (Ridge curvature $\kappa$ dot plot & red mean ± SD)**: Digitized at 1200 dpi (`ry = 1.20 px`):
  - `fig8o_curvature_mean_sd` (`d_fig8o_curvature_mean_sd`, `error_bars: sd`): `1=control` **`0.01692 ± 0.00984 µm⁻¹`**, `2=cyto-D` **`0.05746 ± 0.04934 µm⁻¹`**.
  - `fig8o_curvature_control_scales` (`d_fig8o_curvature_control_scales`): all **12/12** individual control dots (`0.00185..0.0369 µm⁻¹`, sample mean `0.0170`, SD `0.0095`).
  - `fig8o_curvature_cytod_scales` (`d_fig8o_curvature_cytod_scales`): all **12/12** individual cyto-D dots (`0.00436..0.1636 µm⁻¹`, sample mean `0.0574`, SD `0.0492`). Exact 12 vs 12 dot count settles `q_fig8_NO_swapped_refs_and_n` (caption misprints `n = 20`), and dividing the Fig. 8O scale-level SDs by $\sqrt{12}$ yields `SE = 0.00284` and `0.01424 µm⁻¹` -> **Welch $t = -2.79$, $\text{df} = 11.87$ (~12), $p = 0.0165$**, directly reproducing printed $t = -2.78, \text{df} = 12, p < 0.05$ and settling `q_cytod_curvature_se_welch_mismatch` (see `q_digitize_fig8_NO_welch_settlement`).

## Pages 12–22 (printed pp. 12–15 + SI pp. 1–7)
- No figures (Methods, References, Author Information, SI Title, Supplementary Results, and Supplementary Methods).

## Page 23 (SI p. 8; caption continued on Page 24) — Supplementary Figure 1: Synchrotron X-ray nanotomography
- **SI Fig. 1A–D**: Left out — 3D reconstructed X-ray nanotomography renderings and YZ cross-sectional slices of one blue and one black *H. sara* scale (scale bars only).
- **SI Fig. 1E ($G_z(y)$ blue ridge correlation curve)**: Digitized on a uniform `40 nm` grid (`0..1640 nm`, 42 points) at 1200 dpi (`rx = 0.35 px`, `ry = 2.33 px`) as `figs1e_corr_blue_ridge` (`d_figs1e_corr_blue_ridge`). Running `correlation_peak_ridge_spacing` (`run_corr_peak_blue_ridge`) returns **`909.57 nm`** vs printed `910.36 nm` (`-0.79 nm`, `-0.09%`; `cmp_corr_peak_blue_ridge` reports `not reproduced` under strict 2-decimal precision `±0.005 nm` because curve-port scalar outputs do not propagate lag discretization).
- **SI Fig. 1F ($G_z(y)$ black ridge correlation curve)**: Digitized on `0..1640 nm` step `40 nm` (42 points) at 1200 dpi (`rx = 0.24 px`, `ry = 2.67 px`) as `figs1f_corr_black_ridge` (`d_figs1f_corr_black_ridge`). Running `correlation_peak_ridge_spacing` (`run_corr_peak_black_ridge`) returns **`949.66 nm`** vs printed `948.37 nm` (`+1.29 nm`, `+0.14%`; `cmp_corr_peak_black_ridge`).
- **SI Fig. 1G ($G_x(z)$ blue crossrib correlation curve)**: Digitized on `0..1640 nm` step `40 nm` (42 points) at 1200 dpi (`rx = 0.35 px`, `ry = 0.00 px`) as `figs1g_corr_blue_crossrib` (`d_figs1g_corr_blue_crossrib`). Because the two lower evenly spaced y-tick labels both round to `'0.999'` below `'1.000'`, the bottom tick is calibrated honestly as `0.99875 ± 0.00025` (propagated into `lo/hi`). Running `correlation_peak_crossrib_spacing` (`run_corr_peak_blue_crossrib`) returns **`484.40 nm`** vs printed `483.76 nm` (`+0.64 nm`, `+0.13%`; `cmp_corr_peak_blue_crossrib`).
- **SI Fig. 1H ($G_x(z)$ black crossrib correlation curve)**: Digitized on `0..1640 nm` step `40 nm` (42 points) at 1200 dpi (`rx = 0.24 px`, `ry = 0.00 px`) as `figs1h_corr_black_crossrib` (`d_figs1h_corr_black_crossrib`). Running `correlation_peak_crossrib_spacing` (`run_corr_peak_black_crossrib`) returns **`605.62 nm`** vs printed `606.81 nm` (`-1.19 nm`, `-0.20%`; `cmp_corr_peak_black_crossrib`).

## Page 24 (SI p. 9) — Supplementary Figure 2
- **SI Fig. 2**: Left out — static frame of a 3D X-ray nanotomography movie of a false-colored blue iridescent scale (no plot axes).

## Page 25 (SI p. 10) — Supplementary Figure 3
- **SI Fig. 3**: Left out — static frame of a 3D X-ray nanotomography movie of a false-colored black non-iridescent scale (no plot axes).

## Page 26 (SI p. 11; caption continued on Page 27) — Supplementary Figure 4: SEM morphological analysis pipeline
- **SI Fig. 4A–C**: Left out — SEM micrographs illustrating length/width transects (A), PeakFinder green peak overlay (B), and FIJI Steger Ridge Detection overlays (C).
- **SI Fig. 4D**: Left out — 2D `(x, y)` pixel-coordinate mask showing the spatial positions of all 11,361 individual ridge-width measurement pixels across the ~20–30 detected ridge tracks (axis ticks are pixel positions without physical units; not a functional curve, and the per-ridge width summaries are digitized in Panel E).
- **SI Fig. 4E (Per-ridge width mean ± SD across Ridges 1–30)**: Digitized all **30/30** black mean circles and vertical SD bars at 1200 dpi (`ry = 1.80 px`, `error_bars: sd`) as `figs4e_ridge_width` (`d_figs4e_ridge_width`), spanning mean ridge widths from `0.2811 µm` (Ridge 20) to `0.3577 µm` (Ridge 18), grand mean `0.3143 µm`.

## Page 27 (SI p. 12)
- No new figures (SI Fig. 4 caption continuation).

## Page 28 (SI p. 13; caption continued on Page 29) — Supplementary Figure 5: Full adult *H. sara* SEM beeswarm & boxplots
- **SI Fig. 5A (Ridge spacing across 4 scale groups, `n = 800`)**: Digitized all 4 medians and full whisker/dot ranges at 1200 dpi (`ry = 0.30 px`, `error_bars: range`) as `figs5a_ridge_spacing_box` (`d_figs5a_ridge_spacing_box`):
  - `1=black cover`: median `0.9536 µm` (IQR `[0.906, 1.022]`, range `[0.783..1.442 µm]`).
  - `2=black ground`: median `0.9447 µm` (IQR `[0.897, 1.001]`, range `[0.764..1.269 µm]`).
  - `3=blue cover`: median **`0.7483 µm`** (IQR `[0.697, 0.840]`, range `[0.622..1.412 µm]`) — closely matching `actin_templated_ridge_spacing` (`29.5 / 40 = 0.7375 µm`) and settling `q_actin_vs_adult_ridge_spacing_offset`!
  - `4=blue ground`: median **`0.8073 µm`** (IQR `[0.722, 0.892]`, range `[0.631..1.278 µm]`).
- **SI Fig. 5B (Scale area across 4 scale groups, `n = 800`)**: Digitized all 4 medians and full ranges at 1200 dpi (`ry = 1.70 px`, `error_bars: range`) as `figs5b_scale_area_box` (`d_figs5b_scale_area_box`): `1=black cover` `3018 [2295..4580]`, `2=black ground` `4107 [2738..5180]`, `3=blue cover` `2694 [1925..3615]`, `4=blue ground` `3689 [2548..5551] µm²` (all matching printed means `3044, 4123, 2700, 3708 µm²`).
- **SI Fig. 5C (Ridge width across 4 scale groups, plotted `n = 400`)**: Digitized all 4 medians and full ranges at 1200 dpi (`ry = 1.50 px`, `error_bars: range`) as `figs5c_ridge_width_box` (`d_figs5c_ridge_width_box`): `1=black cover` `0.3053 [0.261..0.347]`, `2=black ground` `0.3029 [0.258..0.333]`, `3=blue cover` `0.3108 [0.267..0.387]`, `4=blue ground` `0.3174 [0.251..0.389] µm`.

## Page 29 (SI p. 14) — Supplementary Figure 6: Adult *Heliconius erato demophoon* non-iridescent control
- **SI Fig. 6A, Ai, B**: Left out — dorsal forewing photographs of *H. e. demophoon* vs *H. sara* and SEM micrograph of *H. e. demophoon* ridges.
- **SI Fig. 6C (Ridge number vs scale width regressions & 95% CI ribbons)**: Digitized both OLS regression lines and shaded 95% CI ribbons at 1200 dpi (`rx = 0.00 px`, `ry = 2.10 px`, `error_bars: ci95`):
  - `figs6c_reg_proximal` (`d_figs6c_reg_proximal`, 9 points `W = 26..62 µm`, line residual `0.99 px`): rises from `29.38` to `68.11`. Sweeping `ridge_number_from_width_and_spacing` (`N = W / 0.906 µm`, `run_figs6c_ridge_num_proximal_sweep`) gives `cmp_figs6c_ridge_num_proximal`: **`reproduced`**.
  - `figs6c_reg_distal` (`d_figs6c_reg_distal`, 8 points `W = 26..60 µm`, line residual `0.80 px`): rises from `28.66` to `60.16`. Sweeping `ridge_number_from_width_and_spacing` (`N = W / 0.956 µm`, `run_figs6c_ridge_num_distal_sweep`) gives `cmp_figs6c_ridge_num_distal`: **`reproduced`**.
- **SI Fig. 6D (Ridge number vs ridge spacing, all 20 individual-by-region-by-type means)**: Digitized all **20/20** completely isolated markers at 1200 dpi (`rx = 1.10 px`, `ry = 1.05 px`):
  - `figs6d_proximal_cover` (`d_figs6d_proximal_cover`, 5 blue circles): spacing `0.8496..0.9405 µm`, ridge number `27.04..43.44`.
  - `figs6d_proximal_ground` (`d_figs6d_proximal_ground`, 5 blue triangles): spacing `0.8711..0.9775 µm`, ridge number `49.84..69.39`.
  - `figs6d_distal_cover` (`d_figs6d_distal_cover`, 5 grey circles): spacing `0.8574..1.0007 µm`, ridge number `25.66..34.00`.
  - `figs6d_distal_ground` (`d_figs6d_distal_ground`, 5 grey triangles): spacing `0.9101..1.0871 µm`, ridge number `46.02..59.01`.
  - Overall 10-point regional means (`0.9056 µm` proximal, `0.9559 µm` distal) match printed `0.906` and `0.956 µm` to `0.0004 µm`, while between-individual spread directly confirms `q_hedem_paired_t_se_pooling_mismatch` (see `q_digitize_figS6D_hedem_se_settlement`).

## Page 30 (SI p. 15) — Supplementary Figure 7
- **SI Fig. 7**: Left out — static frame of an animated TauSTED Z-stack of actin bundles at 44% development (no plot axes).

## Page 31 (SI p. 16) — Supplementary Figure 8: Additional TauSTED images of developing *H. sara* scales
- **SI Fig. 8A–F, Ai–Fi**: Left out — merged and depth-color-coded TauSTED fluorescence micrographs at 56%, 69%, and 81% development (no plot axes).

## Page 32 (SI p. 17; caption ends on Page 33) — Supplementary Figure 9: TauSTED actin–chitin colocalization in *P. arcas*
- **SI Fig. 9A–E**: Left out — merged TauSTED micrographs, maximum intensity projection, and YZ orthogonal views (no plot axes).
- **SI Fig. 9F (YZ intensity profile across splayed *P. arcas* ridges, `0.00–12.90 µm`)**: Digitized both curves on a `0.15 µm` grid (87 points each) at 1200 dpi (white tick marks on black background, `rx = 1.83 px`, `ry = 1.90 px`):
  - `figs9f_actin` (`d_figs9f_actin`, green phalloidin): four inter-wall peaks centered near `1.1`, `3.6–4.2`, `7.2–7.5`, and `10.5 µm` (heights `~900–1300`).
  - `figs9f_chitin` (`d_figs9f_chitin`, magenta CBD): four double-walled chitin ridge pairs peaking up to `~4100–4600`, bracketing each actin peak.

## Page 33 (SI p. 18)
- No new figures (SI Fig. 9 caption end only).

## Page 34 (SI p. 19; caption on Page 35) — Supplementary Figure 10: Individual wing reflectance spectra (printed rotated 90° CW)
- **SI Fig. 10A (6 × 6 individual wavelength spectra, `350–700 nm`, rotated 90° CCW upright)**: Calibrated all 6 facet rows independently on y-ticks `[25, 20, 15, 10, 5]%` at 900 dpi (`ry = 0.30–0.85 px`):
  - `figs10a_control_peaks` (`d_figs10a_control_peaks`): all **15/15** red `GIM_TREAT_*` control facets (`11,657–14,598` red px, `0` teal px); peaks `5.77%..16.55%`, mean **`11.24 ± 3.87%`**.
  - `figs10a_treated_peaks` (`d_figs10a_treated_peaks`): all **21/21** teal `TREATED_*` cyto-D facets (`0` red px, `10,907–14,000` teal px — including `TREATED_1B, 1E, 1F`); peaks `2.59%..23.90%` (single outlier `TREATED_6B` at `23.90%`), mean **`7.86 ± 4.65%`**.
  - Settles `q_audit_figS10_orientation_reading` over `q_figS10_panel_coloring_mismatch`, and shows baseline-subtracted peak amplitude `max - min` (`7.30 ± 2.77%` vs `2.97 ± 3.41%`, Welch $t = 4.19, \text{df} = 33.35, p = 0.00019 < 0.001$) closely reproduces printed $t = 4.34, \text{df} = 33, p < 0.001$ (clarifying `q_pavo_reflectance_metric_unstated`; see `q_digitize_figS10_color_and_pavo_metric`).
- **SI Fig. 10B (6 × 6 individual angular reflectance dots, `0°–20°`, rotated 90° CCW upright)**: Calibrated all 6 facet rows independently on y-ticks `[25, 20, 15, 10, 5, 0]%` at 900 dpi (`ry = 0.36–0.73 px`):
  - `figs10b_control_angle_peaks` (`d_figs10b_control_angle_peaks`): all **15/15** red `GIM_TREAT_*` facets (`18,446–18,539` red px, `0` teal px); peaks `5.75%..16.57%`, mean **`11.23 ± 3.88%`** (matching Panel A individual-by-individual to `±0.02%`!).
  - `figs10b_treated_angle_peaks` (`d_figs10b_treated_angle_peaks`): all **21/21** teal `TREATED_*` facets (`0` red px, `18,782–19,032` teal px); peaks `2.11%..23.91%`, mean **`7.54 ± 4.68%`**.

## Page 35 (SI p. 20) — Supplementary Figure 11: Additional SEMs of cytochalasin D treated scales
- **SI Fig. 11A–F**: Left out — SEM micrographs of scale pinching, filled windows, uneven ridge layering, wobbling ridges, and toppled ridges (scale bars only).

## Page 36 (SI p. 21) — Supplementary Figure 12: AFM imaging and ridge height analysis
- **SI Fig. 12A (Treated scale 3D AFM & transverse height profile)**: 3D rendering and 2D inset left out (no plot axes). Transverse height profile digitized for **Profile 4 (blue curve)** on `0.04..2.88 µm` step `0.04 µm` (72 points) at 1800 dpi (`rx = 1.82 px`, `ry = 2.50 px`) as `figs12a_afm_treated_profile4` (`d_figs12a_afm_treated_profile4`); Profiles 1–3 overlap heavily in dark/green/red so the cleanest unoccluded profile spanning the full 3 µm scan was read. Running `afm_ridge_height_from_profile` (`run_afm_height_treated_p4`) returns **`H_mean = 0.7553 µm`**, **`H_max = 0.9000 µm`**.
- **SI Fig. 12B (Control scale 3D AFM & transverse height profile)**: 3D rendering and 2D inset left out. Transverse height profile digitized for **Profile 4 (blue curve)** on `0.05..4.50 µm` step `0.05 µm` (90 points) at 1800 dpi (`rx = 2.78 px`, `ry = 2.50 px`) as `figs12b_afm_control_profile4` (`d_figs12b_afm_control_profile4`). Running `afm_ridge_height_from_profile` (`run_afm_height_control_p4`) returns **`H_mean = 1.1006 µm`**, **`H_max = 1.3167 µm`**.
- **SI Fig. 12C (AFM ridge height dot plot & red mean ± SD, mislabeled `µm` on axis)**: Digitized at 1200 dpi (`ry = 0.33 px`) preserving the **literal printed unit `um`** (`q_figS12C_unit_mislabel_um_vs_nm`, `q_digitize_figS12_afm_unit_settlement`):
  - `figs12c_ridge_height_mean_sd` (`d_figs12c_ridge_height_mean_sd`, `error_bars: sd`): `1=control` **`1304.2 ± 165.0`**, `2=cyto-D` **`1121.3 ± 187.2`** (literally labeled `µm`, numerically `nm`).
  - `figs12c_ridge_height_control_dots` (`d_figs12c_ridge_height_control_dots`): all **6/6** control dots (`1092, 1201, 1201, 1373, 1444, 1516`; sample mean `1304.5`, SD `164.5`).
  - `figs12c_ridge_height_cyto_dots` (`d_figs12c_ridge_height_cyto_dots`): all **6/6** cyto-D dots (`906.8, 1001, 1001, 1144, 1277, 1395`; sample mean `1120.8`, SD `187.1`). Direct comparison with `run_afm_height_control_p4` (`1.317 µm`) and `run_afm_height_treated_p4` (`0.900 µm`) confirms the 1000× `µm`-for-`nm` mislabeling.

## Page 37 (SI p. 22)
- No figures (Supplementary References continuation).
