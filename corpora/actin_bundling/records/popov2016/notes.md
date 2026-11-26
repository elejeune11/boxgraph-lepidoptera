# Cross-Record `relate` & `validate` Audit — `popov2016` (MEDYAN)

Owner record: `popov2016` (Popov, Komianos & Papoian, *PLoS Comput. Biol.* 12(4):e1004877, 2016).
Audited pairs (4 remaining ready pairs; `chandrasekaran2019+popov2016` and `maxian2022+popov2016` are handled in those partner records):
1. `bathe2008 + popov2016` (WLB vs MEDYAN)
2. `chandrasekaran2024 + popov2016` (CytoSim / LAMMPS vs MEDYAN)
3. `freedman2019 + popov2016` (AFINES vs MEDYAN)
4. `graham2023 + popov2016` (Discrete-Ring VASP Droplet vs MEDYAN)

---

## Summary Counts (inside `popov2016/`)

- **`cited.yaml`**: **0** links (`popov2016` cites zero corpus records; `bathe2008` is the only older paper among the 4 pairs and is not cited in `popov2016`'s 125 main + SI references).
  - *Partner-side citation flagged (read-only partner, not edited here)*: `graham2023` SI Table 1 (p. 39) cites `popov2016` as `ref-si-10` for `single_filament_k_stretch = 100 pN/nm` (corresponding to `popov2016/K_actin_str = 100 pN/nm`, Table 2B p. 60), and `freedman2019` cites `popov2016` qualitatively as `ref-37` (p. 5) for alternative simulation frameworks. `graham2023/references.yaml` (`ref-si-10`) currently lacks `source: popov2016` and `graham2023/cited.yaml` lacks the `single_filament_k_stretch -> popov2016/K_actin_str` link.
- **`relations.yaml`**: **15** cross-record relations added (**0** `bathe2008`, **7** `chandrasekaran2024`, **5** `freedman2019`, **3** `graham2023`).
- **`validations.yaml`**: **0** entries (`popov2016` has 0 `origin: measured` datums, and none of the 62 `origin: measured` datums across the 4 partners matches a runnable output quantity/dimension of `popov2016`).
- **`questions.yaml`**: 4 cross-framework & audit questions appended (`q_cross_framework_stretching_and_wlb_coupling`, `q_cross_framework_brownian_ratchet_and_slip_bond`, `q_cross_framework_validations_audit`, `q_filament_bending_float64_cancellation_note`).

---

## Pair-by-Pair Audit

### 1. `bathe2008 + popov2016` (WLB vs MEDYAN)
- **Citations / Relations / Validations**: 0 / 0 / 0.
- **Why 0 relations**: All 31 `bathe2008` models compute continuum/FE Wormlike Bundle flexural rigidity ($\kappa_B$ in `N*m^2`), inter-filament coupling $\alpha$, shear modulus $G_B$, Euler buckling force, fracture density, or Langmuir spacing. `popov2016`'s 49 models describe single-filament discrete potentials, mechanochemical rate laws, and reaction-diffusion parameter conversions; bundles emerge only in the non-runnable 3D simulator `medyan_3d_actomyosin_simulator`.
- **Why 0 validations**: All 6 `origin: measured` datums in `bathe2008` are Claessens 2006 / Fig. 4 bundle bending stiffnesses $\kappa_B$ (`N*m^2`), which no `popov2016` model predicts.
- **Sharp WLB vs MEDYAN comparison**:
  - Both share Kojima et al. (1994) (`bathe2008/EfAf_actin = 44 nN`, `popov2016/ref-98`) and single-filament persistence length (`kappa_f ≈ 7.0e-26 N*m^2` in `bathe2008` vs `l_p k_B T = 7.26e-26 N*m^2` in `popov2016`).
  - MEDYAN deliberately softens `K_actin_str = 100 pN/nm`, giving effective extensional rigidity $E_f A_f^{\text{MEDYAN}} = K_{\text{str}} l_{\text{cyl}} = 10.8\,\text{nN}$ ($l_{\text{cyl}}=108\,\text{nm}$) or $2.7\,\text{nN}$ ($l_{\text{cyl}}=27\,\text{nm}$) — **4.1×–16.3× softer** than physical actin. Inside Bathe's WLB coupling parameter $\alpha = k_x L^2 / (E_f A_f \delta)$, this shifts effective $\alpha$ upward by **4.1×–16.3×** toward the fully-coupled regime.
  - WLB assumes permanent crosslinkers with $k_\perp \sim 1000\,\text{pN/nm}$ and shear stiffness $k_x \sim 10^{-5}\text{--}10^{-6}\,\text{N/m}$, whereas MEDYAN treats $\alpha$-actinin as central-force Hookean springs ($K_\alpha = 8\,\text{pN/nm}$, 125× softer than $k_\perp$, zero intrinsic angular/shear stiffness) that stochastically unbind via a Bell slip bond ($k_0 = 0.3\,\text{s}^{-1}$, $F_u = 17.2\,\text{pN}$).

### 2. `chandrasekaran2024 + popov2016` (CytoSim / LAMMPS vs MEDYAN)
- **Citations / Relations / Validations**: 0 / 7 / 0.
- **Relations verified (7)**:
  - `rel_chandrasekaran2024_bell_slip_bond` (`same_as`, symbolic): `popov2016/alpha_actinin_slip_bond_s3eq6` ↔ `chandrasekaran2024/bell_slip_bond_unbinding` (`output_map: {k_unbind: k_unbind}`).
  - `rel_chandrasekaran2024_brownian_ratchet_scales` (`scales_as`): `popov2016/actin_brownian_ratchet_s3eq7` → `chandrasekaran2024/brownian_ratchet_growth_rate` (identical exponential Brownian ratchet up to monomer increment $l_{\text{mon}} = 2.7\,\text{nm}$).
  - `rel_chandrasekaran2024_lammps_bond_scales_filament_stretch` (`scales_as`): `popov2016/filament_stretching_energy_eq2` → `chandrasekaran2024/lammps_harmonic_bond_potential`.
  - `rel_chandrasekaran2024_lammps_bond_scales_crosslinker_stretch` (`scales_as`): `popov2016/motor_crosslinker_stretching_energy_eq9` → `chandrasekaran2024/lammps_harmonic_bond_potential`.
  - `rel_chandrasekaran2024_lammps_angle_scales_filament_bend` (`scales_as`, small angle $\theta \in [0.01, 0.2]\,\text{rad}$): `popov2016/filament_bending_energy_eq1` → `chandrasekaran2024/lammps_harmonic_angle_potential`.
  - `rel_chandrasekaran2024_kuramoto_scales_brownian_noise` (`scales_as`): `popov2016/kuramoto_length_s2eq1` → `chandrasekaran2024/brownian_noise_stddev` (3D $\sqrt{6 D \tau}$ vs 1D $\sqrt{2 D dt}$, exact dimensionality factor $\sqrt{3}$).
  - `rel_chandrasekaran2024_molecule_count_same_as_nmiia_count` (`same_as`, numeric `< 1e-15`): `popov2016/nmiia_minifilament_count_from_ratio` ↔ `chandrasekaran2024/molecule_count_from_concentration` (cubic volume vs equivalent-volume sphere, CODATA $N_A$).
- **Why 0 validations**: All 10 `origin: measured` datums in `chandrasekaran2024` are VASP-droplet disk/ring/shell percentages, absent from MEDYAN.
- **Sharp CytoSim/LAMMPS vs MEDYAN comparison**:
  - **6.67× Brownian ratchet stall-force divergence**: `popov2016` derives $F_{\text{actin,poly}} = k_B T / l_{\text{mon}} = 1.5\,\text{pN}$ (Peskin 1993 / Footer 2007), whereas `chandrasekaran2024` sets $f_{\text{grow}} = 10\,\text{pN}$ — at $5\,\text{pN}$ load, MEDYAN drops growth to $3.6\%$ vs CytoSim's $60.7\%$ (**17× divergence**).
  - **Shared Bell slip bond vs catch bond**: Both MEDYAN $\alpha$-actinin ($F_u = 17.2\,\text{pN}$) and CytoSim VASP ($f_u = 10\,\text{pN}$) obey Bell's slip bond, opposite in load sign to MEDYAN's NMIIA catch bond (`nmiia_unbinding_and_attachment_time_s3eq2`).
  - **Inverse-variable F-actin conservation**: `popov2016/mean_filament_length_mass_balance` solves for $L_{\text{mean}}$, the algebraic inverse of `chandrasekaran2024/final_actin_concentration` solving for $C_{\text{actin}}$.

### 3. `freedman2019 + popov2016` (AFINES vs MEDYAN)
- **Citations / Relations / Validations**: 0 (in `popov2016`; `freedman2019` cites `popov2016` qualitatively as `ref-37`, p. 5) / 5 / 0.
- **Relations verified (5)**:
  - `rel_freedman2019_filament_stretching_same_as` (`same_as`, symbolic): `freedman2019/filament_potential_si_eq1` (`U_stretch`) ↔ `popov2016/filament_stretching_energy_eq2` (`U_str`).
  - `rel_freedman2019_filament_bending_derived_small_angle` (`derived_from`, symbolic `series: theta -> 0, 3`): `popov2016/filament_bending_energy_eq1` → `freedman2019/filament_potential_si_eq1` (`U_bend`).
  - `rel_freedman2019_filament_bending_approx_medyan` (`approximated_by`, max err `0.53%` over $\theta \in [0.01, 0.35]\,\text{rad}$): `freedman2019/filament_potential_si_eq1` (`U_bend`) ≈ `popov2016/filament_bending_energy_eq1`.
  - `rel_freedman2019_crosslinker_stretching_same_as` (`same_as`, numeric `1.64e-15`): `freedman2019/crosslinker_potential_si_eq4` (`U_stretch`) ↔ `popov2016/motor_crosslinker_stretching_energy_eq9` (`U_motor`).
  - `rel_freedman2019_gap_bending_approx_two_medyan_hinges` (`approximated_by`, max err `0.24%` over $\Delta l_{xl} \in [5, 35]\,\text{nm}$, $l_g \in [100, 1000]\,\text{nm}$): `freedman2019/gap_bending_energy_eq1` (`U`) ≈ two MEDYAN cosine hinges (`2 * eps_bend`).
- **Why 0 validations**: 6 of `freedman2019`'s 8 measured datums are crosslinker domain lengths (`um`); the other 2 are TIRF actin polymerization speeds (`22.15 ± 3.78 nm/s` at `0.75 uM` and `46.99 ± 3.62 nm/s` at `1.5 uM`), whereas `popov2016` has no runnable model outputting actin polymerization velocity in `[length]/[time]`. (Notably, multiplying `popov2016` Table 2A's Fujiwara net barbed-end rate `11.6*[A] - 1.4 s^-1` by `l_mon = 2.7 nm` yields `19.71 nm/s` and `43.20 nm/s` — within `0.65 SD` and `1.05 SD` of Freedman's TIRF measurements.)
- **Sharp AFINES vs MEDYAN comparison**:
  - **Time-integration vs spring stiffness trade-off**: AFINES uses explicit Brownian dynamics ($\Delta t = 2 \times 10^{-5}\,\text{s}$), forcing ultra-soft springs (`k_a = 5 pN/um` at `l_a = 1 um` → $E_f A_f = 5\,\text{pN}$, **2,160× softer** than MEDYAN at `108 nm` and **8,800× softer** than Kojima 1994; `k_xl = 10 pN/um`, **800× softer** than MEDYAN's `8 pN/nm`). MEDYAN decouples chemistry from quasi-static CG minimization (`g_tol = 1 pN` every `0.05 s`), permitting stiff springs.
  - **Crosslinker thermodynamics**: AFINES enforces equilibrium detailed balance via Metropolis MC on `U_xl = U_stretch + U_bend + U_bind` (with perpendicular bending stiffness $\kappa_{xl}^B$), whereas MEDYAN uses non-equilibrium Bell slip-bond unbinding on pure central-force springs.

### 4. `graham2023 + popov2016` (Discrete-Ring VASP Droplet vs MEDYAN)
- **Citations / Relations / Validations**: 0 (in `popov2016`; `graham2023` cites `popov2016` as `ref-si-10` for `single_filament_k_stretch = 100 pN/nm`) / 3 / 0.
- **Relations verified (3)**:
  - `rel_graham2023_discrete_bending_constant_derived` (`derived_from`, symbolic `subs: N_fil = 1`): `graham2023/discrete_bending_constant` → `popov2016/actin_cylinder_bending_energy_table2b`.
  - `rel_graham2023_ring_confinement_derived_medyan_bend` (`derived_from`, symbolic `series: theta -> 0, 3`): `popov2016/filament_bending_energy_eq1` → `graham2023/ring_confinement_bending_energy`.
  - `rel_graham2023_ring_confinement_approx_medyan_bend` (`approximated_by`, max err `0.21%` over $L \in [0.05, 0.15]\,\mu\text{m}$, $R \in [0.5, 4.0]\,\mu\text{m}$): `graham2023/ring_confinement_bending_energy` ≈ `popov2016/filament_bending_energy_eq1`.
- **Why 0 validations**: All 38 `origin: measured` datums in `graham2023` measure VASP droplet FRAP, diameters, aspect ratios, partition coefficients, morphology abundances, or droplet-confined bundle elongation (`23 ± 3 nm/s` at `2 uM` actin), none matching a runnable `popov2016` model.
- **Sharp Graham vs MEDYAN comparison**:
  - Graham et al. directly adopt MEDYAN's cylindrical-segment discretization (`L_cyl = 100 nm` vs `108 nm`, cosine hinge `k_bend*(1 - u_i . u_{i+1})` with `k_bend = N_fil * L_p * k_B T / L_cyl`, and `k_stretch = N_fil * 100 pN/nm` citing Popov 2016 `[10]`), but drop the `1/2` stretching prefactor in SI Eq. (S2) (`sum k_stretch*(||r_i - r_{i+1}|| - L_cyl)^2`), effectively doubling per-filament stretching stiffness to `200 pN/nm` ($E_f A_f = 20\,\text{nN}$, still 2.2× softer than physical actin).

---

## Cross-Framework Parameter Comparison Matrix

| Physical Parameter | Physical / WLB (`bathe2008`) | MEDYAN (`popov2016`) | Graham Ring (`graham2023`) | AFINES (`freedman2019`) | CytoSim / LAMMPS (`chandrasekaran2024`) |
|---|---|---|---|---|---|
| Single-filament $l_p$ / $\kappa_f$ | $\kappa_f \approx 0.070\,\text{pN}\cdot\mu\text{m}^2$ ($17.0\,\mu\text{m}$) | $l_p = 17.7\,\mu\text{m}$ ($0.0726\,\text{pN}\cdot\mu\text{m}^2$) | $L_p = 17.73\,\mu\text{m}$ ($0.0729\,\text{pN}\cdot\mu\text{m}^2$) | $L_p = 17\,\mu\text{m}$ (KMC); $\kappa_B = 0.01\text{--}0.5$ (AFINES) | $\kappa = 0.075\,\text{pN}\cdot\mu\text{m}^2$ ($18.25\,\mu\text{m}$) |
| Segment length $l_{\text{seg}}$ | Continuum (crosslink spacing $\delta$) | $27\,\text{nm}$ / $108\,\text{nm}$ | $100\,\text{nm}$ | $1000\,\text{nm}$ ($1\,\mu\text{m}$) | CytoSim constraint; LAMMPS $0.1\sigma$ |
| Filament spring $K_{\text{str}}$ | $k_f = E_f A_f/\delta$ ($407\,\text{pN/nm}$ at $108\,\text{nm}$) | $100\,\text{pN/nm}$ (with $\tfrac{1}{2}$) | $100\,\text{pN/nm}$ (no $\tfrac{1}{2}$ → $200\,\text{pN/nm}$ eff.) | $0.005\,\text{pN/nm}$ ($5\,\text{pN/\mu m}$, with $\tfrac{1}{2}$) | Inextensible (CytoSim); $30000\,k_B T/\sigma^2$ (LAMMPS) |
| Effective $E_f A_f = K_{\text{eff}} l_{\text{seg}}$ | $44,000\,\text{pN}$ (Kojima 1994) | $10,800\,\text{pN}$ ($108\,\text{nm}$) / $2,700\,\text{pN}$ ($27\,\text{nm}$) | $20,000\,\text{pN}$ | $5\,\text{pN}$ ($8,800\times$ softer) | $\infty$ (CytoSim) |
| Discrete bending potential | Continuum $\tfrac{1}{2}\kappa_B \int \kappa^2 ds$ | $\varepsilon_{\text{bend}}(1 - \cos\theta)$ | $k_{\text{bend}}(1 - \mathbf{u}_i\cdot\mathbf{u}_{i+1})$ | $\tfrac{\kappa_B}{2 l_a}\theta^2$ (small-angle) | $K_{\text{angle}}(\theta - \pi)^2$ (LAMMPS, $2500\,k_B T$) |
| Crosslinker stiffness | $k_\perp \sim 1000\,\text{pN/nm}$; $k_x \sim 10^{-5}\text{--}10^{-6}\,\text{N/m}$ | $K_\alpha = 8\,\text{pN/nm}$ (central force) | Implicit (decoupled $k_b \propto N_{\text{fil}}$) | $k_{xl} = 0.01\,\text{pN/nm}$ + $\kappa_{xl}^B$ bending | Harmonic spring (CytoSim/LAMMPS) |
| Crosslinker unbinding | Permanent (fixed) | Bell slip bond ($F_u = 17.2\,\text{pN}$, $k_0 = 0.3\,\text{s}^{-1}$) | Permanent (continuum bundle) | Detailed-balance Metropolis MC ($k_{\text{off}} = 0.05\,\text{s}^{-1}$) | Bell slip bond ($f_u = 10\,\text{pN}$, $k_0 = 0.1\text{--}1.75\,\text{s}^{-1}$) |
| Brownian ratchet force | N/A | $F_{\text{actin,poly}} = 1.5\,\text{pN}$ ($k_B T / l_{\text{mon}}$) | N/A | Force-free Markov growth ($k_g \Delta l_a$) | $f_{\text{grow}} = 10\,\text{pN}$ ($6.67\times$ larger) |
