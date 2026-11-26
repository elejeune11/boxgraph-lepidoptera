# Figure Inspection — `kierfeld2005`

Direct inspection of every page and figure panel of Kierfeld, Kühne & Lipowsky, *Phys. Rev. Lett.* **95**, 038102 (2005) (`corpus/actin_bundling/2005_Kierfeld.pdf`, pages 1–4). All four pages were rendered at 300 dpi with pixel rulers (`figwork/p1.png`–`figwork/p4.png`) via `boxgraph figure page`, and both quantitative plot panels in Fig. 3 were rendered at **1200 dpi** (`figwork/fig3a_1200.png`, `figwork/fig3b_1200.png`).

**Summary:** 14 calibrated `mode: code` readings were authored in `figures.yaml` covering all 7 branches of **Fig. 3(a)** and all 7 branches of **Fig. 3(b)**, and every reading was visually verified against its magenta overlay (`figwork/fig3a_*.png`, `figwork/fig3b_*.png`). **Fig. 1** (3D Monte Carlo snapshots) and **Fig. 2** (schematic diagram) carry no numeric axes or scales and are left out.

---

## Page-by-page inspection

- **Page 1 (`038102-1`, `figwork/p1.png`) — Fig. 1(a, b, c) — Left out:**
  - *Visual contents:* Three color 3D Monte Carlo snapshots of $N = 20$ filament bundles: (a) loose bundle close to the unbinding transition dominated by pair contacts, (b) kinetically trapped segregated sub-bundles deep in the bound phase, (c) equilibrium roughly cylindrical hexagonal bundle.
  - *Reason left out:* Qualitative 3D snapshots with no coordinate axes, scale bar, or numeric values.

- **Page 2 (`038102-2`, `figwork/p2.png`) — Fig. 2 — Left out:**
  - *Visual contents:* Schematic drawing of two parallel filaments displaced along $\mathbf{z}_1(x)$ and $\mathbf{z}_2(x)$ along axis $x$, bridged by flexible cross-linkers with two adhesive end groups.
  - *Reason left out:* Purely schematic illustration; no numeric axes or plotted data.

- **Page 3 (`038102-3`, `figwork/p3.png`, `figwork/fig3a_1200.png`, `figwork/fig3b_1200.png`) — Fig. 3(a) and Fig. 3(b) — Digitized (14 calibrated readings):**
  - **Fig. 3(a)** (`crop: [1700, 8140, 5080, 10700]` at `1200 dpi`, calibration residuals `x_res = 0.36 px`, `y_res = 0.34 px`):
    - Plots mean total energy per filament per length $\langle \mathcal{H} \rangle / (N L)$ (in units of $T$, with $\Delta x = 1$) vs. effective potential strength $|\bar{W}| \in [4, 9]$ (note minor EPS font quirk where decimal points in half-integer tick labels `4.5, 5.5, 6.5, 7.5, 8.5` are not visible).
    - All 7 branches digitized and overlay-verified:
      - `fig3a_N2` ($N=2$, diamonds, 14 points spanning unbound plateau $+1.0$ and bound branch $|\bar{W}| \in [5.8, 8.0]$; vertical jump at $|\bar{W}_c| \approx 5.70$).
      - `fig3a_N3` ($N=3$, pluses, 14 points; vertical jump at $|\bar{W}_c| \approx 5.31$).
      - `fig3a_N5` ($N=5$, open squares, 13 points; vertical jump at $|\bar{W}_c| \approx 5.01$).
      - `fig3a_N10_compact` ($N=10$, compact cylindrical lower branch, asterisks, 10 points; vertical jump at $|\bar{W}_c| \approx 4.70$; markers end at $|\bar{W}| \approx 6.2$, with only a dotted guide continuing to $8.0$).
      - `fig3a_N20_compact` ($N=20$, compact cylindrical lower branch, open diamonds, 14 points; vertical jump at $|\bar{W}_c| \approx 4.61$; descends to $-12.51$ at $|\bar{W}| = 8.0$).
      - `fig3a_N20_planar` ($N=20$, planar-initial segregated upper branch, triangles with thick solid line, 9 points on $|\bar{W}| \in [4.8, 8.0]$).
      - `fig3a_N10_planar` ($N=10$, planar-initial segregated upper branch, crosses with thick solid line, 7 points on $|\bar{W}| \in [5.5, 8.0]$).
  - **Fig. 3(b)** (`crop: [5140, 8140, 8560, 10700]` at `1200 dpi`, log-log axes calibrated from frame borders + major `0.1` ticks, residuals `x_res = 0.33 px`, `y_res = 0.17 px`):
    - Plots mean interfilament separation $\langle \Delta z \rangle \equiv \langle |\Delta \mathbf{z}_{ij}| - \ell_r \rangle$ (in units of $\Delta x$) vs. reduced potential strength $(|\bar{W}| - |\bar{W}_c|)/|\bar{W}| \in [0.01, 1]$.
    - All 7 branches digitized and overlay-verified: `fig3b_N2` (9 points), `fig3b_N3` (8 points), `fig3b_N5` (12 points), `fig3b_N10_compact` (9 points), `fig3b_N20_compact` (10 points), `fig3b_N10_planar` (9 points), `fig3b_N20_planar` (10 points).
    - *Direct visual findings from Fig. 3(b):* (1) In the equilibrium compact cylindrical branch deep in the bound phase ($x \gtrsim 0.1$), $\langle \Delta z \rangle$ monotonically **increases** with $N$ ($N=2 < N=3 < N=5 < N=10 < N=20$, from $\sim 0.012$ to $\sim 0.19$), directly contradicting the text on p. 4 col 1 stating that "the bundle thickness and the mean energy per filament of the equilibrium form decrease with increasing $N$, as can be seen in Fig. 3" (only the mean energy per filament decreases with $N$, and $\langle \Delta z \rangle$ decreases with $|\bar{W}|$, not with $N$). (2) The digitized $N=2$ values ($\sim 0.0118$ at $|\bar{W}|=8.0$) are ~5–6× larger than $\langle |\Delta\mathbf{z}_{12}| - \ell_r \rangle \approx 0.0019$ from equilibrium Monte Carlo of the stated Hamiltonian (which matches Fig. 3(a) energy to ~1%), and instead closely match the RMS separation $\sqrt{\langle (|\Delta\mathbf{z}_{12}| - \ell_r)^2 \rangle} \approx 0.011$ (logged in `questions.yaml`).

- **Page 4 (`038102-4`, `figwork/p4.png`) — No figures:**
  - Continuation of MC discussion, Discussion & Conclusion, Acknowledgments, and References [1]–[17].

---

## Dedicated Digitization Verification Pass (`digitize:kierfeld2005`, 2026-11-25)

- **Completeness across all 4 pages (`pages/page-001.png`–`page-004.png`):** Independently inspected all four pages. Confirmed **Fig. 1(a–c)** (p. 1, qualitative 3D MC snapshots without axes or scale bars) and **Fig. 2** (p. 2, schematic illustration of two cross-linked filaments) carry no quantitative data; **Page 4** has no figures. Confirmed **Fig. 3(a)** and **Fig. 3(b)** (p. 3) are the sole quantitative panels, each containing 7 MC data branches (`N = 2, 3, 5`, `N = 10` compact/planar, `N = 20` compact/planar), with no error bars and no insets. All 14 branches are digitized (**148 points total**: 81 across Fig. 3(a) and 67 across Fig. 3(b)); no printed `~` values in `data.yaml` are eye-read figure values.
- **Script speed & calibration check:** Verified all 14 `figures/*.py` readers operate purely on the passed numpy `image` array (zero PDF re-rendering; 0.47–1.72 s per `read(image, at)` call on the 1200-dpi crops; 148/148 points found, 0 errors). Calibration residuals verified from the frame and interior ticks: `x_res = 0.36 px, y_res = 0.34 px` (Fig. 3(a), linear axes) and `x_res = 0.33 px, y_res = 0.17 px` (Fig. 3(b), log10 axes).
- **Overlay regeneration & visual inspection:** Re-ran `boxgraph figure check` on all 14 readings (refreshing overlays including the post-audit fixes to `fig3a_N5.py` and `fig3a_N10_compact.py`), assembled four compact 2×2 contact sheets (`figwork/contact_fig3a_1.png`, `figwork/contact_fig3a_2.png`, `figwork/contact_fig3b_1.png`, `figwork/contact_fig3b_2.png`, width 1400 px) plus high-zoom verification sheets (`figwork/zoom_fig3a.png`, `figwork/zoom_fig3b.png`, `figwork/zoom_remaining.png`, `figwork/zoom_check_3b.png`), and visually verified every picked mark. All 81 Fig. 3(a) marks sit cleanly on their target markers/thick solid lines (with `fig3a_N10_compact` stopping at `|W_bar| = 6.1` where markers end and only a dotted guide continues). All 67 Fig. 3(b) marks sit cleanly on their target markers/curves, with conservative overlap-widened intervals where curves cross or crowd (`fig3b_N2` at `x >= 0.10`, `fig3b_N3` at `x = 0.19, 0.24`, and `fig3b_N20_compact` at `x = 0.23` where `N = 10` planar crosses `N = 20` compact).
- **Saved readings:** Generated and saved `readings.json` (all 14 readings, 148 points) via `boxgraph.engine.save_readings` and validated `SOUND` with 0 errors.

