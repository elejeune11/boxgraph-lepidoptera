# Figure Inventory and Digitization Log — `totz2024`

Complete panel-by-panel audit of every page of `corpus/actin_bundling/2024_Totz.pdf` (23 PDF pages: pp. 1–14 main article & Experimental Procedures, pp. 15–23 Supplemental Information). All 47 quantitative readings were implemented as 900–1200 dpi `mode: code` scripts under `figures/`, run through `boxgraph figure check`, and visually inspected via their magenta overlay PNGs in `figwork/`.

## Page 1 — Graphical Abstract
- **Graphical Abstract (p. 1)**: Left out — composite journal summary graphic re-rendering thumbnails of Figs. 1A, 1E–G, 2C, 3A/C/E, and 3H; introduces no new quantitative data beyond the main figures.

## Page 2 (printed p. 1)
- No figures (Title, Abstract, Introduction text only).

## Page 3 (printed p. 2) — Figure 1: *In vivo* phase data of developing *Vanessa cardui* scales
- **Fig. 1A (i–iv)**: Left out — butterfly photograph (i), optical scale micrograph (ii), and SEM images of mature cover scales and ridges (iii–iv) with scale bars only; no quantitative axes.
- **Fig. 1B–D**: Left out — 2D quantitative phase surface-height heatmaps (`-50` to `+150 nm`) at 35.8% (B), 40.0% (C), and 41.0% (D) pupal development; their quantitative 1D centerline profiles are digitized directly in Figs. 1E–G.
- **Fig. 1E (35.8% centerline profile)**: Digitized as `fig1E_profile_35_8pct` (`d_fig1E_profile_35_8pct`, 1200 dpi, `res_x = 0.78 px`, `res_y = 0.00 px`, 19 points `x = 0.25–4.75 um`; smooth surface fluctuating between ~0 and 12 nm).
- **Fig. 1F (40.0% centerline profile)**: Digitized as `fig1F_profile_40_0pct` (`d_fig1F_profile_40_0pct`, 1200 dpi, `res_x = 0.90 px`, `res_y = 0.00 px`, 19 points `x = 0.25–4.75 um`; nascent protoridges peaking at ~25–68 nm).
- **Fig. 1G (41.0% centerline profile)**: Digitized as `fig1G_profile_41_0pct` (`d_fig1G_profile_41_0pct`, 1200 dpi, `res_x = 0.90 px`, `res_y = 0.67 px`, 19 points `x = 0.25–4.75 um`; regular protoridges peaking at ~68–95 nm).

## Page 4 (printed p. 3)
- No figures (Main text only).

## Page 5 (printed p. 4) — Figure 2: Constrained-growth plate buckling model
- **Fig. 2A**: Left out — fluorescence projection micrograph of parallel actin bundles and scale membrane adapted from Pomerantz et al. (2021) (microscopy image, no plotted axes).
- **Fig. 2B**: Left out — TEM cross-section of actin bundles and developing cuticle/membrane adapted from Pomerantz et al. (2021) (microscopy image, no plotted axes).
- **Fig. 2C**: Left out — qualitative schematic diagram of growth, compression, bending, actin-membrane clamping, and cytoplasmic pressure.
- **Fig. 2D**: Left out — uncalibrated 3D numerical plate-buckling rendering without numeric coordinate axes.
- **Fig. 2E — Growth strain $g$ subpanel (top-left)**: Digitized all 3 curves at 1200 dpi (`res_x = 0.00 px`, `res_y = 0.60 px`, 7 points each at `x = -0.4..+0.4 um`): `fig2E_g_1_0pct` (`d_fig2E_g_1_0pct`, peak `208.9 nm`), `fig2E_g_0_5pct` (`d_fig2E_g_0_5pct`, peak `160.3 nm`), and `fig2E_g_0_0pct` (`d_fig2E_g_0_0pct`, peak `102.9 nm`). Confirms all profiles plot `w(x) = a*(1 + cos(2*pi*x/L))` (peak `2*a = 208.5, 160.1, 102.4 nm` under corrected Eq. 3, omitting the factor of 1/2).
- **Fig. 2E — Young's modulus $Y$ subpanel (top-right)**: Digitized all 3 curves at 1200 dpi (`res_x = 1.21 px`, `res_y = 0.60 px`, 7 points each): `fig2E_Y_1e8` (`d_fig2E_Y_1e8`, peak `133.2 nm`), `fig2E_Y_1e6` (`d_fig2E_Y_1e6`, peak `160.3 nm`), and `fig2E_Y_1e5` (`d_fig2E_Y_1e5`, peak `259.1 nm`). Matches corrected Eq. 3 `2*a` at `g = 0.5%` (`132.9, 160.1, 259.6 nm`), not the caption's `g = 1.0%`.
- **Fig. 2E — Pressure $p$ subpanel (bottom-right)**: Digitized all 3 curves at 1200 dpi (`res_x = 1.40 px`, `res_y = 0.80 px`, 7 points each): `fig2E_p_100` (`d_fig2E_p_100`, peak `259.7 nm`), `fig2E_p_50` (`d_fig2E_p_50`, peak `217.4 nm`), and `fig2E_p_10` (`d_fig2E_p_10`, peak `159.9 nm`). Matches corrected Eq. 3 `2*a` at `g = 0.5%` (`259.6, 217.5, 160.1 nm`), not the caption's `g = 1.0%`.
- **Fig. 2E — Bundle distance $L$ subpanel (bottom-left)**: Digitized all 3 curves at 1200 dpi (`res_x = 0.00 px`, `res_y = 0.80 px`): `fig2E_L_1_5` (`d_fig2E_L_1_5`, 7 points, peak `127.9 nm`), `fig2E_L_1_0` (`d_fig2E_L_1_0`, 7 points, peak `74.4 nm`), and `fig2E_L_0_5` (`d_fig2E_L_0_5`, 5 points, peak `11.3 nm`). Matches corrected Eq. 3 `2*a` at `g = 0.5%` with base spacing `L_0 = 1.0 um` (`128.1, 74.8, 10.6 nm`), not the caption's `g = 1.0%` and `L = 1.8 um`.

## Page 6 (printed p. 5)
- No figures (Main text and Eqs. 3–4 only).

## Page 7 (printed p. 6) — Figure 3: Validation of mechanical buckling hypothesis
- **Fig. 3A (39.6% development)**: Digitized yellow dashed model profile (`fig3A_model_39_6pct` / `d_fig3A_model_39_6pct`, 1200 dpi, `res_x = 1.00 px`, `res_y = 0.60 px`, 9 dashes, least-squares cosine amplitude `w_0 = 24.31 nm`) and *in vivo* black mean curve with gray ±1 SD band (`fig3A_invivo_39_6pct` / `d_fig3A_invivo_39_6pct`, `error_bars: sd`, 9 unobstructed segments).
- **Fig. 3B (40.0% development)**: Digitized orange dashed model profile (`fig3B_model_40_0pct` / `d_fig3B_model_40_0pct`, 1200 dpi, `res_x = 1.40 px`, `res_y = 0.60 px`, 10 dashes, `w_0 = 33.58 nm`) and *in vivo* mean ± 1 SD (`fig3B_invivo_40_0pct` / `d_fig3B_invivo_40_0pct`, `error_bars: sd`, 9 segments).
- **Fig. 3C (40.3% development)**: Digitized red dashed model profile (`fig3C_model_40_3pct` / `d_fig3C_model_40_3pct`, 1200 dpi, `res_x = 0.80 px`, `res_y = 0.60 px`, 11 dashes, `w_0 = 41.69 nm`) and *in vivo* mean ± 1 SD (`fig3C_invivo_40_3pct` / `d_fig3C_invivo_40_3pct`, `error_bars: sd`, 10 segments).
- **Fig. 3D (40.7% development)**: Digitized magenta dashed model profile (`fig3D_model_40_7pct` / `d_fig3D_model_40_7pct`, 1200 dpi, `res_x = 0.60 px`, `res_y = 0.60 px`, 11 dashes, `w_0 = 42.31 nm`) and *in vivo* mean ± 1 SD (`fig3D_invivo_40_7pct` / `d_fig3D_invivo_40_7pct`, `error_bars: sd`, 9 segments).
- **Fig. 3E (41.0% development)**: Digitized purple dashed model profile (`fig3E_model_41_0pct` / `d_fig3E_model_41_0pct`, 1200 dpi, `res_x = 0.60 px`, `res_y = 0.60 px`, 14 dashes, `w_0 = 64.31 nm`) and *in vivo* mean ± 1 SD (`fig3E_invivo_41_0pct` / `d_fig3E_invivo_41_0pct`, `error_bars: sd`, 13 segments).
- **Fig. 3F (scaled amplitude $w_0$ violin plots)**: Digitized means (`fig3F_w0_mean` / `d_fig3F_w0_mean`, 1200 dpi, `res_y = 1.20 px`: `31.90, 36.68, 53.41, 48.46, 72.08 nm`) and medians with Q1–Q3 interquartile range bars (`fig3F_w0_median_iqr` / `d_fig3F_w0_median_iqr`, `error_bars: range`: medians `31.05, 39.29, 54.50, 52.92, 71.93 nm`). Translucent swarm scatter of 3,722 individual profiles left out (overplotted raw points summarized by mean, median, and IQR).
- **Fig. 3G (unscaled width $L$ violin plots)**: Digitized means (`fig3G_L_mean` / `d_fig3G_L_mean`, 1200 dpi, `res_y = 0.00 px`: `1.277, 1.595, 1.501, 1.440, 1.598 um`) and medians with Q1–Q3 IQR bars (`fig3G_L_median_iqr` / `d_fig3G_L_median_iqr`, `error_bars: range`: medians `1.282, 1.642, 1.467, 1.519, 1.566 um`). Raw swarm scatter left out.
- **Fig. 3H (phase-space amplitude isocontours $g(Y)$)**: Digitized all 5 dashed amplitude isocontours at 900 dpi (`res_x = 0.20 px`, `res_y = 0.20 px`): `fig3H_iso_39_6pct` (7 dashes, `g = 0` intercept at `Y = 3.70 x 10^6 Pa`), `fig3H_iso_40_0pct` (9 dashes), `fig3H_iso_40_3pct` (7 dashes), `fig3H_iso_40_7pct` (7 dashes), and `fig3H_iso_41_0pct` (16 dashes). Continuous background heatmap encodes the same analytical Eq. 3 surface and is left out in favor of the 5 vector dashed isocontours.

## Page 8 (printed p. 7)
- No figures (Main text only).

## Page 9 (printed p. 8) — Figure 4: Inter-bundle spacing across butterfly species
- **Fig. 4A**: Left out — TEM cross-section of developing *Colias eurytheme* scale adapted from Ghiradella (1974) with wing photograph inset from Ficarrotta et al. (2022) (micrograph/photograph; derived spacings $L_{1,A}, L_{2,A}$ are read in Fig. 4C).
- **Fig. 4B**: Left out — confocal fluorescence micrograph of *Agraulis vanillae* silver scale adapted from Dinwiddie et al. (2014) with wing inset (micrograph/photograph; derived spacings $L_{1,B}, L_{2,B}$ are read in Fig. 4C).
- **Fig. 4C (phase space of amplitude vs bundle spacing $L$ and growth strain $g$)**: Digitized at 1200 dpi (`res_x = 1.40 px`, `res_y = 0.60 px`, calibrating `x` from interior major ticks `0.5, 1.0, 1.5 um` and right spine `2.0 um` because the left spine sits at `L = 0.01 um`):
  - `fig4C_species_spacings` (`d_fig4C_species_spacings`): all 5 vertical black dashed species lines (`L_2,A = 0.3406 um`, `L_2,B = 0.3599 um`, `L_1,A = 0.8006 um`, `L_1,B = 0.9098 um`, `L_Vcardui = 1.601 um`).
  - `fig4C_iso_10nm` (`d_fig4C_iso_10nm`, 5 points along `L = 0.51–0.75 um`), `fig4C_iso_1nm` (`d_fig4C_iso_1nm`, 5 points along `g = 0.25–0.65%`), and `fig4C_iso_0_1nm` (`d_fig4C_iso_0_1nm`, 5 points along `g = 0.05–0.68%`): all 3 gray amplitude isocontours. Continuous background heatmap left out in favor of the exact isocontours.

## Pages 10–14 (printed pp. 9–13)
- No figures (Discussion, Experimental Procedures, Acknowledgments, Author Contributions, Declaration of Interests, References).

## Page 15 (Supplemental Information cover)
- No figures.

## Page 16 (caption continued on p. 17) — Figure S1: *In vivo* phase unwrapping and segmentation pipeline
- **Fig. S1A**: Left out — 2D wrapped phase microscopy image (`[-pi, pi]`), no quantitative plot axes.
- **Fig. S1B**: Left out — 2D compass-colormap phase-gradient orientation image for manual ROI selection.
- **Fig. S1C**: Left out — 2D cropped ROI wrapped phase image.
- **Fig. S1D**: Left out — 3D perspective rendering of unwrapped ROI height (5x height exaggeration); quantitative cross-sections are plotted directly below in S1E.
- **Fig. S1E (multi-cross-sectional ROI view)**: Digitized at 900 dpi (`res_x = 0.67 px`, `res_y = 0.20 px`): red mean height profile (`figS1E_mean_profile` / `d_figS1E_mean_profile`, 11 unobstructed points `x = 0.5–9.5 um`), 5 blue peak markers (`figS1E_peaks` / `d_figS1E_peaks`), and 6 green trough markers (`figS1E_troughs` / `d_figS1E_troughs`, bounding 5 ridges of widths `1.42, 1.42, 1.86, 1.42, 1.77 um`). Underlying translucent gray individual cross-section spaghetti left out.
- **Fig. S1F (5 individually segmented and rotated protoridges)**: Digitized the red mean height profile of all 5 subplots at 1200 dpi (`res_x = 0.20–0.40 px`, `res_y = 0.40 px`): `figS1F_ridge1` (`d_figS1F_ridge1`), `figS1F_ridge2` (`d_figS1F_ridge2`), `figS1F_ridge3` (`d_figS1F_ridge3`), `figS1F_ridge4` (`d_figS1F_ridge4`), and `figS1F_ridge5` (`d_figS1F_ridge5`). Underlying gray individual trace spaghetti left out.
- **Fig. S1G**: Left out — 5 width-normalized subplots containing only overlapping translucent gray individual traces with no mean or summary curve.
- **Fig. S1H**: Left out — single-ROI collection of overlapping translucent gray normalized ridge traces with no summary curve.
- **Fig. S1I**: Left out — conceptual stacked thumbnail graphic with ellipsis (`...`) illustrating pooling across ROIs for one developmental time point.
- **Fig. S1J (rescaled all-ridge catalog at one time point)**: Digitized red mean curve and blue ±1 SD curves at 1200 dpi (`figS1J_mean_sd` / `d_figS1J_mean_sd`, `error_bars: sd`, `res_x = 0.50 px`, `res_y = 0.15 px`, 7 points `x = 0.2–1.4 um`; peak `61.75 nm, SD [28.84, 94.56] nm`, matching 41.0% development in Fig. 3E). Underlying gray individual trace cloud left out.

## Pages 17–23 (Supplemental Table S1, Notes S1–S2, Supplemental References)
- No figures (Table S1 on p. 17 is already transcribed in `data.yaml`; pp. 17–23 contain Supplemental Notes S1–S2 equations and references only).
