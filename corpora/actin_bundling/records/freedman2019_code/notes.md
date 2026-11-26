# Audit Notes: `freedman2019_code` (`Simfreed/AFINES@growing_spacer_exv`)

- **Repository**: `https://github.com/Simfreed/AFINES/tree/growing_spacer_exv`, commit `0c6c349`
- **Tree digest**: `sha256:bd1357a511addd53b246a09d02a430cc047cafac865b035fb53b276c5455b03b`
- **Companion paper**: `freedman2019` — Freedman et al., *PNAS* **116**(33), 16192–16197 (2019)
- **Additional corpus link**: `winkelman2016` — Winkelman et al., *Curr. Biol.* **26**(20), 2697–2706 (2016) (embedded `expt_data_txt` table in `xlink_sorting_paper/kmc.ipynb` Cell 9)

---

## 1. Repository Architecture & Runnable Inventory

| Component | Source files | BoxGraph model(s) | Status |
|---|---|---|---|
| 1D KMC bending & growing-bundle model | `xlink_sorting_paper/kmc.ipynb` (Cells 0–12) | `kmc_notebook_bending_and_switching`, `kmc_notebook_growing_bundle` | Live `software: true` adapters executing unmodified notebook code via `$BOXGRAPH_SOURCE` |
| AFINES C++ parameter & diffusivity mapping | `prog/bundles.cpp`, `src/motor_ensemble.cpp:313-317`, `src/bead.cpp:25`, `src/spacer.cpp:55`, `xlink_sorting_paper/config/{grow,nogrow}.cfg`, `README.md`, `process/plots.nb` | `afines_parameter_mapping` | Exact Python port |
| AFINES C++ barbed-end growth & splitting | `src/filament.cpp:673-736` | `afines_barbed_end_growth_cpp` | Exact IEEE-754 Python port |
| AFINES C++ link-link excluded volume | `src/filament_ensemble.cpp:571-715`, `src/spring.cpp:261-358` | `afines_excluded_volume_cpp` | Exact Python port |
| AFINES C++ filament mechanics & BD step | `src/spring.cpp:77-81,376-378`, `src/filament.cpp:80-83,206-242,518-575` | `afines_filament_mechanics_cpp` | Exact Python port |
| AFINES C++ spacer mechanics & Metropolis attachment | `src/spacer.cpp:208-388`, `src/motor.cpp:486-491`, `src/motor_ensemble.cpp:313-317` | `afines_spacer_mechanics_and_mc_cpp` | Exact Python port (both buggy C++ and physical energy change) |
| Mathematica lattice & domain post-processing | `process/afines2lattice.m`, `process/domain_length.m`, `process/cytomod_functions.m:858-914`, `process/plots.nb:767-847` | `afines_lattice_and_domain_processing` | Exact Python port |
| Integrated 2-filament AFINES BD trajectory | `prog/bundles.cpp:410-485` + C++ force/integrator routines | `afines_two_filament_bundle_port` | Short-horizon (0.01 s) stochastic Python port, `dt`-refined |
| Full Boost C++ `bin/bun` 2000–6000 s sweep + missing `lattice.py` | `prog/bundles.cpp`, `xlink_sorting_paper/README.md` | `afines_bundles_cpp_driver` | `runnable: false` (Boost build, hundreds of CPU-h, missing `lattice.py` & `crashed_list.txt`) |

**Record totals**: 59 quantities, 98 data (61 cross-linked to `freedman2019`), 9 added probes, 10 models, 17 runs, 29 comparisons (**15 reproduced, 2 close, 7 not reproduced, 5 no verdict**), 2 settled refinements (`ref_kmc_trials`, `ref_bundle_dt`, both `belongs_to: model`), 11 open questions, 8 references.

---

## 2. Five Core Code-vs-Paper Discoveries

### 2.1 Fig. S3A KMC 39% Discrepancy 100% Solved (`xlink_sorting_paper/kmc.ipynb`)
- **Growth algorithm**: Cell 5 (`extension_population_ratio`) implements a **continuous-time Gillespie (Poisson) process**, NOT deterministic periodic site addition. `site_extension_rate = float(extension_rate)/site_width` is appended as the final rate in `tmp_rate_list` alongside binding rates `kon1`, `p_diff*kon1` (or `p_diff*kon2`, `kon2`) over all `n_sites_away` vacant sites, and time advances by `-np.log(r2)/tmp_cum_rate[-1]`.
- **Root cause of the 39% discrepancy**: Cell 2 (lines 98–100) **transposes 35 nm (alpha-actinin size) and 37 nm (F-actin helical half-pitch)**:
  - `site_width = 35.` (yielding `max_n_sites = int(15000 / 35) = 428`, vs the paper's stated `l_b = 37 nm` and `N = 405`).
  - `dxl = 37. - 8. = 29.` nm (vs the paper's stated `Delta l_xl = 35 - 8 = 27 nm`).
  - Because `f(x, Lp) = (Lp*1e3/site_width)/x * arcsin((dxl/site_width)/x)^2` scales approximately as `dxl^2 / site_width^3`, swapping `(37, 27) -> (35, 29)` inflates the bending barrier by `(29/27)^2 * (37/35)^3 = 1.362` (+36.2%) at every site.
- **Exact reproduction**: Executing unmodified `kmc.ipynb` (`run_kmc_figs3a_code_exact`) with Cell 0's `np.random.seed(1)` + Cell 5's trailing 1-trial warmup + Cell 6's `num_trials = 1000` at `(site_width = 35 nm, dxl = 29 nm)` yields:
  - `N_diff / N_same = 0.016472` (`%.3f` -> **`0.016`** printed in Fig. S3A),
  - `epsilon = 4.1061` (`%.1f` -> **`4.1`** printed in Fig. S3A).
- **Paper parameters check**: Running the **exact same unmodified function** (`run_kmc_figs3a_paper_params`) at the paper's printed `(37 nm, 27 nm)` yields `0.022624` and `3.7887` (+41% higher than `0.016`, matching `freedman2019`'s independent Gillespie check `0.0222` / `3.81`).
- **Equal-binding gap**: Cell 4's `20` sites × `35 nm` = `0.700 um` (`U = 0.0417 k_B T`, `p_switch = 0.4896`), explaining the paper's `~0.7 um` (whereas `20 * 37 nm = 0.740 um`).
- **Additional `kmc.ipynb` bugs**:
  1. Cells 9–12 (lines 797, 869, 910, 942, 950) pass `std_domain_lengths[0]` / `std_domain_lengths_mu[0]` (Type A / fascin error bars) for **both** fascin and alpha-actinin `plt.errorbar` calls.
  2. Cell 5 computes domain lengths from instantaneous filament tip length (`filament_length`) at swap events rather than bound cross-linker site position (`current_site * site_width`).
  3. Cell 10 line 876 packages Cell 9's `avg_domain_lengths` into `competition_prediction_data` instead of `avg_domain_lengths_mu`.
- **Cell 9 `kratio` sweep (`cmp_kmc_cell9_eps_curve`)**: Running `kmc_notebook_growing_bundle` at `num_trials = 150` × 3 replicates reproduces 4/6 points within 2 SE (`kratio = 1.0, 1.5, 2.0, 2.5`) and agrees within 0.2%–1.8% on `kratio = 0.3` (`4.426` vs `4.417`) and `0.5` (`4.010` vs `4.082`), which flag `not reproduced` because Cell 9's raw 15-digit stdout floats carry no printed MC SE (~0.04 at 1000 trials) and the 3-replicate sample SE happened to be tiny (`0.0017` and `0.0050`).

### 2.2 SI Eq. 3 Excluded-Volume Units & Crossing Bugs (`src/filament_ensemble.cpp:571-715`)
- **Units resolved**: Line 673 computes `pe_exv += kexv * pow(1 - r*b, 2)` (`b = 1/rmax`) and lines 668–671 compute `Fx1 = 2*kexv*dist[0]*b*(1/r - b)`, which equals exact `-grad U` (verified in `run_exv_fig2d`: `U = 0.02 pN*um`, `F_code = F_analytic = 1.6 pN` at `kexv = 0.08`, `rmax = 0.05`, `r = 0.025`). Therefore **`kexv` has energy units `pN*um`** (`0.08 pN*um = 20 k_B T` at `temperature = 0.004`), NOT `pN*um^2` as misprinted in Table S1 & Fig. 2D.
- **Config/README omission**: Both `grow.cfg:37` and `nogrow.cfg:37` set `kexv = 0`, and `README.md` omits `--kexv 0.08` for Fig. 2D, even though `plots.nb:8457,8460` proves Fig. 2D ran with `kexv = 0.08`.
- **Midpoint-crossing zero-force bug**: The outer guard `if (r < rmax)` (line 618) checks only endpoint-to-segment distances `r_c[0..3]`. Two 1-um links crossing orthogonally at their midpoints have `r = 0.5 um > 0.05 um`, skip the entire repulsion block, and experience **`0.0 pN`** excluded-volume force (`exv_F_cross_mid = 0.0`).
- **Near-end crossing 45° doubled-force bug**: When intersecting links do satisfy `r < rmax` (lines 700–713), C++ applies a hardcoded +45° force `Fx1 = Fy1 = 2*kexv*b/sqrt(2)` to **both** endpoints of link 1 and `(-Fx1, -Fy1)` to **both** endpoints of link 2, ignoring segment orientation and doubling total repulsive force to `4*kexv/rmax = 6.4 pN` (`exv_F_cross_end = 6.4`).

### 2.3 SI Section A Bending Calibration, 5× `set_bending` Division & the `(200, 300) nm` Gap
- **Undocumented 5× division (`src/motor_ensemble.cpp:313-317`)**: `spacer_ensemble::set_bending(modulus, ang)` divides the config parameter by spacer rest length: `double kb = modulus / mld;`. Thus `s1_bend = 0.005` (`spacer1_length = 0.2 um`) yields an actual angular stiffness `kb_short = 0.025 pN*um` in `U_xl^bend = 0.5*kb*(alpha - pi/2)^2` — **5× larger** than Table S1's printed `0.005 pN*um` (`cmp_kb_short_eff_vs_table_s1` not reproduced).
- **Four mutually inconsistent long cross-linker bending prescriptions** (leaving the `(200, 300) nm` gap, including `227 nm`, unresolved):
  1. **SI Section A**: `kappa_long^B / kappa_short^B = 22*(l_long/l_short) - 21` for `l_long <= 200 nm`, and `4` for `l_long >= 300 nm` (undefined on `(200, 300) nm`).
  2. **`grow.cfg:28` / `nogrow.cfg:28`**: `s2_bend = 0.06` at `l_long = 0.3 um` -> raw ratio `0.06/0.005 = 12` (the exact value of `22*(300/200) - 21 = 12`), effective ratio `kb_long / kb_short = (0.06/0.3)/(0.005/0.2) = 0.20/0.025 = 8.0` (both `not reproduced` vs SI's `4`).
  3. **`xlink_sorting_paper/README.md:40,102` (Fig. 2C)**: `kb0 = 0.04`, `kb2 = kb0 * l2` -> constant effective `kb_long = 0.04 pN*um` (effective ratio `1.6`), passed via invalid flag `--s2bend` (`prog/bundles.cpp:135` requires `--s2_bend`).
  4. **`process/plots.nb:14229` (Fig. 5C, 5E, S5, Movies S6–S7)**: loads job directories labeled `_s2kb0.0675_` -> effective `kb_long = 0.0675/0.3 = 0.225 pN*um` (effective ratio `9.0`, raw ratio `13.5`).

### 2.4 SI Eq. 2 Barbed-End Split Condition & IEEE-754 Remainder (`src/filament.cpp:673-736`)
- C++ tests `if ( lb + dL < l0_max )` — splitting at `>= 2*l_a`, NOT the strict `> 2*l_a` printed below SI Eq. (2) — and resets both split links to `spring_l0 = 1.0 um`.
- Exact divisors `{0.01, 0.02, 0.04, 0.05, 0.10} um` discard zero remainder (`100.0 nm/s` at `0.10 um`, `10` steps/split, `t_done = 140 s`).
- Non-divisors in Fig. 5C–E still lose remainder: `lgrow = 0.03 um` discards `0.02 um`/split (`34` steps -> `29.41 nm/s`, `-2.0%`), and `lgrow = 0.07 um` discards `0.05 um`/split (`15` steps -> `66.67 nm/s`, `-4.8%`).
- IEEE-754 binary double accumulation causes `1.0 + 400*0.0025 = 1.9999999999999998 < 2.0`, adding 1 extra step at `0.0025 um` (`401` steps -> `2.4938 nm/s`, `t_done = 5614 s`) and `0.005 um` (`201` steps -> `4.9751 nm/s`).
- Commit `0c6c349` ("updated default l0max") changed `prog/network.cpp:189` to `default_value(2)` but left `prog/bundles.cpp:178` at `default_value(0)` (overridden by `l0max=2` in `grow.cfg`/`nogrow.cfg`).

### 2.5 Table S1 vs Code Protocol Mismatches, C++ Bugs & Hidden Defaults
- **Initial length (`grow.cfg:7`)**: `nmonomer = 2` (= **1 link = `1.0 um`**, confirmed by `plots.nb:14255`'s `tDoneGrowing = Ceiling[14/lgrow]`), matching Fig. 5B caption (`1 um`) and showing Table S1's `N_{l,0} = 2, 15` confused bead count (`2`) with link count (`15` from `nogrow.cfg`'s `nmonomer = 16`).
- **Minimum growth rate (`README.md:256`, `plots.nb:14195`)**: Fig. 5C starts at `lgrow = 0.0025 um` (`2.5 nm/s`, matching Movie S6), proving Table S1's `0.008 um` (`8 nm/s`) is a misprint.
- **Unstated `tfinal` & `dt` overrides (`plots.nb`)**: Slow growth (`0.0025, 0.005`) ran to `tfinal = 6000 s` (`plots.nb:14217`, 3× Table S1's `2000 s`); Fig. 5B `s2len = 0.3` used `dt = 0.00001 s` (`plots.nb:11266,11270`, ½ Table S1's `0.00002 s`).
- **Hardcoded `temperature = 0.004 pN*um`**: Used in `grow.cfg`, `nogrow.cfg`, `prog/bundles.cpp:98`, `README.md:201`, and `plots.nb:8466` instead of `k_B*(300 K) = 0.0041419 pN*um`, giving exact `L_p = {2.5, 17.0, 125.0} um` in Movie S3 and `mu = -0.008/0.004 = -2 k_B T`.
- **C++ `s2_koff` typo (`prog/bundles.cpp:372,376`)**: Passes `s2_kend, s2_kend` instead of `s2_koff, s2_kend` to `small_xlinks` (silently ignoring `s2_koff` whenever `s2_koff != s2_kend`).
- **C++ `spacer::metropolis_prob()` detailed-balance bugs (`src/spacer.cpp:220-224,355-378`)**: Uses pre-attachment `(hx[hd], hy[hd])`, old `disp`, and stale `pos_a_end[hd] = 0` instead of `newpos`; includes only one head's bending energy upon entering `(1,1)`; and `spacer::update_force()` never resets `b_eng` on unbinding, subtracting stale `b_eng[hd]` from `delE`. At our probe (`run_spacer_mc_short_bug`), `delE_buggy = 0.01225` vs `delE_true = 0.01481 pN*um`, inflating attachment probability by 1.90× (`1.87e-6` vs `9.86e-7`).
- **Hidden Stokes radii & numerical defaults**: Bead Stokes radius `R = actin_length = 0.5 um` (`src/bead.cpp:23-25`, `D_bead = 0.4244 um^2/s`, matching Leimkuhler-Matthews `bd_prefactor`); spacer Stokes radius `R = mld` (`src/spacer.cpp:55,151`, giving `D_short = 1.061` vs `D_long = 0.707 um^2/s`, 1.5× species-dependent diffusivity); capture cutoff `s1_cut = s2_cut = 0.063 um` (`prog/bundles.cpp:120,133`); `grid_factor = 3` in `.cfg` vs `2.5` in Table S1; `quad_update_period = 10`; `maxSmallAngle = 0.001`.
- **Broken `README.md` & missing files**: Missing `lattice.py` and `crashed_list.txt`; wrong filename `config/no_grow.cfg`; invalid CLI flags `--s2bend` and `--s1_koff1`; undefined `${lgrow}` and `${s1dens}`; omitted `--spacer2_length` in Fig. 5A–B and `--kexv 0.08` in Fig. 2D; `y0 = 0.375` outside filament `y = ±0.125`. `plots.nb` also reveals lattice MC sampled `t = 9000..10000`, plotted `Delta l_xl = 0` at `10 nm` on Fig. 2C's log axis, and omitted `Delta l_xl = 0.3 um` ("forgot to do 0.3 for lattice -- FIX").
- **Mathematica post-processing verified consistent with SI Section B**: `afines2lattice.m` uses `l0 = 0.037 um` (`405` bins, `Sign[bcs1-bcs2]/.{-1->2}` maps ties to `0`), `cytomod_functions.m` implements `LatticeDomainMerge(nc=4, rc=1 um)`, and `plots.nb` uses `angThresh = 0.8` and domain length `x_last - x_first`.
