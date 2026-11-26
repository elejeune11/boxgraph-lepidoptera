# Cross-Record Relate & Validate Notes — `chandrasekaran2019`

Summary of ModelAtlas (`relations.yaml`), cited-value linkage (`cited.yaml`), cross-record validation (`validations.yaml`), and parameter/equation drift audits (`questions.yaml`) for owner record `chandrasekaran2019` (Chandrasekaran, Upadhyaya & Papoian, *PLoS Comput. Biol.* 15(7): e1007156, 2019) across its four partner pairs.

## 1. `chandrasekaran2019` ↔ `popov2016` (MEDYAN ↔ MEDYAN)

- **Citation link**: `chandrasekaran2019` directly runs and documents `popov2016`'s MEDYAN force field and mechanochemical cycle (`ref-58` / `ref-S6`).
- **`cited.yaml` (5 links)**:
  - `compartment_size_medyan` (`500 nm`, already `from_source: ref-S6`) → `popov2016/l_comp_table2d` (`500 nm`) — exact match.
  - `l0_cylinder` (`108 nm`) → `popov2016/l_cyl_large` (`108 nm`) — exact match.
  - `alpha_actinin_binding_distance` (`30 to 40 nm`) → `popov2016/alpha_actinin_binding_range` (`35 ± 5 nm`) — exact match.
  - `myosin_mf_binding_distance` (`175 to 225 nm`) → `popov2016/nmiia_binding_range` (`200 ± 25 nm`) — exact match.
  - `myosin_stall_force` (`12.6 pN`) → `popov2016/F_nmiia_stall` (`15 pN`) — exposes table-row transcription error (`12.6 pN` is `F_NMIIA,unbind` on the row directly above `F_stall = 15 pN` in Popov Table 2C, p. 61).
  - *Intentionally unlinked*: `gamma_myosin_bound_heads_printed` (`0.05` dimensionless) and `epsilon_myosin_unbind_printed` (`0.2 /pN`) have transposed units relative to Popov's `nmiia_gamma` (`0.05 /pN`) and `nmiia_beta` (`0.2` dimensionless); primary-source datums (`l_p_actin`, `myosin_duty_ratio`, `myosin_hill_zeta`, `myosin_heads_per_mf_sim`) keep their existing `from_source` in `data.yaml`.
- **`relations.yaml` (15 verified cross-record relations)**:
  - `rel_popov2016_mesoscopic_rate_eqS2_corrected` (`same_as`, proposed): `popov2016/compartment_bimolecular_rate_table2a` ↔ `chandrasekaran2019/mesoscopic_rate_eqS2_corrected` (`V_r = l_comp^3`, CODATA `N_A`).
  - `rel_popov2016_actin_stretching_eqS3` (`same_as`, stated): Eq. (S3) ↔ Popov Eq. (2).
  - `rel_popov2016_actin_bending_eqS4` (`same_as`, stated): Eq. (S4) ↔ Popov Eq. (1).
  - `rel_popov2016_excluded_volume_eqS6` (`same_as`, stated): Eq. (S6) ↔ Popov Eq. (4) (both share dropped Jacobian `|l_i||l_j|`, `K_vol` in `pN*nm^5`).
  - `rel_popov2016_excluded_volume_eqS5_jacobian` (`same_as`, proposed): Eq. (S5) ↔ Popov Eq. (4) via `K_vol -> K_vol * L^2`.
  - `rel_popov2016_minifilament_stretching_eqS7` (`same_as`, stated): Eq. (S7) ↔ Popov Eq. (9).
  - `rel_popov2016_crosslinker_stretching_eqS8` (`same_as`, stated): Eq. (S8) ↔ Popov Eq. (9) (crosslinker potential).
  - `rel_popov2016_boundary_repulsion_eqS9` (`same_as`, stated): Eq. (S9) ↔ Popov Eq. (11).
  - `rel_popov2016_minifilament_binding_eqS10` (`derived_from`, `subs: alpha = 1`, proposed): Popov S3 Eq. (1) → Eq. (S10).
  - `rel_popov2016_myosin_bound_heads_eqS11` (`same_as`, proposed): Popov S3 Eq. (3) ↔ `myosin_bound_heads_eqS11_inverse_force` (`N_0 = rho * N_total`).
  - `rel_popov2016_myosin_zero_load_bound_heads` (`derived_from`, `subs: F_ext = 0`, proposed): Popov S3 Eq. (3) → `myosin_zero_load_bound_heads`.
  - `rel_popov2016_myosin_catch_bond_eqS12` (`same_as`, proposed): `myosin_unbinding_catch_bond_eqS12_dimensionless_eps` ↔ Popov S3 Eq. (2) (`beta = epsilon`).
  - `rel_popov2016_myosin_unloaded_walking_eqS13` (`same_as`, stated): Popov S3 Eq. (5) ↔ Eq. (S13) (`N_bound_0 = rho * N_total`).
  - `rel_popov2016_myosin_hill_walking_eqS14` (`same_as`, proposed): Popov S3 Eq. (4) ↔ Eq. (S14) under `F_stall -> N_total * F_stall` (proving algebraically that Eq. S14 absorbs `N_total` into `F_stall`, making printed `12.6 pN` 17.9×–35.7× below Popov's ensemble stall force `225–450 pN`).
  - `rel_popov2016_alpha_actinin_slip_bond_eqS15` (`same_as`, stated): Eq. (S15) ↔ Popov S3 Eq. (6).
- **Sharp drift & bug propagation logged (`q_popov2016_medyan_parameter_and_equation_drift`)**:
  - **Propagated upstream bugs**: (i) dropped Jacobian `|l_i||l_j|` in Eq. (3)→(4) / Eq. (S5)→(S6) (with `K_vol` printed `10^5 pN/nm^4` in Popov Table 2B vs omitted in Chandrasekaran); (ii) off-by-one nematic sum `i = 0..N` divided by `N` in Popov Eq. (14) → Chandrasekaran Eq. (S18) (both contradicted by S2 Fig `S <= 0.992`).
  - **New regressions in `chandrasekaran2019`**: (i) inverted bimolecular rate conversion `k_mu * (V_r / N_A)` in Eq. (S2); (ii) corrupted Kuramoto length formula on p. 29; (iii) `/pN` unit swapped from `γ = 0.05/pN` onto `ϵ = 0.2` (Popov `β = 0.2`) between Eqs. (S11)–(S12); (iv) `F_stall = 12.6 pN` copied from Popov `F_NMIIA,unbind` while dropping `1/N_total`; (v) omitted `d_step / l_site = 6/27` step-size conversion in Eq. (S13).
  - **Simulation protocol drifts**: `D_alpha/D_actin` increased 2.5× (`1/10 → 1/4`) and `D_NMIIA/D_actin` 12.5× (`1/100 → 1/8`); binding sites changed from dual occupancy to mutual exclusion (4 exclusive sites / 108 nm); minifilament heads changed from `10–30` to `15–30`; motor stoichiometry shifted from `R_m:a = 0.005–0.02` (molecules:actin) to `M:A = 0.0225–0.675` (heads:actin, up to 33.75× higher, unlocking aster collapse at `M:A >= 0.225`). Consistent: `χ = 1.0` treadmilling rates in S11 Fig (`0.46–0.58 monomers/(fil. s)`) match Popov's `0.5 s^-1`; S10 Fig plots normalized `R_g(t)/R_g(0)` avoiding Popov's absolute `R_g = 1.77 um > sqrt(3)*L/2` paradox.
- **Validations**: 0 (both records are pure simulation papers with `0` `origin: measured` datums).

## 2. `chandrasekaran2019` ↔ `chandrasekaran2024` (MEDYAN 2019 vs CytoSim/LAMMPS 2024)

- **`cited.yaml`**: 0 links (2019 precedes 2024).
- **`relations.yaml` (3 verified proposed relations)**:
  - `rel_chandrasekaran2024_bell_slip_bond_eqS15` (`same_as`): exact Bell slip-bond unbinding identity (`crosslinker_slip_bond_unbinding_eqS15` ↔ `chandrasekaran2024/bell_slip_bond_unbinding`, `F_u = 17.2 pN` α-actinin vs `10 pN` VASP tetramer).
  - `rel_chandrasekaran2024_actin_stretching_scales_lammps_bond` (`scales_as`): MEDYAN Eq. (S3) (`pN*nm`) ↔ LAMMPS harmonic bond potential (`k_B T`, `sigma = 30 nm, k_B T = 4.11 pN*nm`).
  - `rel_chandrasekaran2024_actin_bending_scales_lammps_angle` (`scales_as`): MEDYAN Eq. (S4) cosine hinge bending ↔ LAMMPS harmonic angle bending (`theta0 = pi`) in the semiflexible small-angle regime `theta ∈ [0.005, 0.1] rad`.
- **Validations**: 0 (`chandrasekaran2019` has `0` measured datums; all 10 `origin: measured` datums in `chandrasekaran2024` are WT/GABmut VASP liquid-droplet shell/ring/disk abundance percentages without motors or α-actinin). Logged in `q_chandrasekaran2024_medyan_vs_cytosim_comparison_and_no_validation`.

## 3. `bathe2008` ↔ `chandrasekaran2019` (WLB Continuum Theory vs MEDYAN Bundles)

- **`cited.yaml`**: 0 links (neither paper cites the other).
- **`relations.yaml`**: 0 relations (orthogonal scales: macroscopic continuum bundle bending stiffness `kappa_B(N, alpha)` vs microscopic discrete segment potentials, stochastic RDME kinetics, and morphology order parameters).
- **Validations**: 0 (`chandrasekaran2019` has `0` measured datums and no model predicting `bathe2008`'s 6 measured Claessens 2006 fascin bundle bending stiffness curves `kappa_B(N)`).
- **Scale diagnostic (`q_bathe2008_orthogonal_scales_and_wlb_regime`)**: On `chandrasekaran2019`'s initial `N = 30, L = 2 um, d_lat = 35 nm` α-actinin bundle (`delta ≈ 33–1350 nm` from S4 Fig, `chi^2 = 48`), Bathe's shear stiffness `k_x = 10^-6 – 10^-5 N/m` gives `alpha = 10^-4 – 0.034` (decoupled/weakly coupled, `kappa_B/(N kappa_f) = 1.0–2.3`), whereas substituting MEDYAN's axial α-actinin spring `K_alpha = 8 pN/nm = 8×10^-3 N/m` would give `alpha = 0.73–27` (`kappa_B/(N kappa_f) = 29–602`).

## 4. `chandrasekaran2019` ↔ `freedman2019` (MEDYAN vs AFINES)

- **`cited.yaml`**: 0 links (`ref-83`, Freedman bioRxiv 493841, is cited qualitatively on p. 7 for passive crosslinker length-mismatch sorting; 0 numerical values imported).
- **`relations.yaml` (4 verified proposed relations)**:
  - `rel_freedman2019_actin_stretching_eqS3` (`same_as`): `freedman2019/filament_potential_si_eq1` (`U_stretch`) ↔ MEDYAN Eq. (S3).
  - `rel_freedman2019_actin_bending_derived_harmonic` (`derived_from`, `series: theta -> 0, 4`): exact small-angle Taylor expansion of MEDYAN Eq. (S4) cosine hinge bending → AFINES `filament_potential_si_eq1` (`U_bend`) under `eps_bend = kappa_B / l_a`.
  - `rel_freedman2019_actin_bending_approx_harmonic` (`approximated_by`, `theta ∈ [0.01, 0.35] rad`, tol `2%`): quantitative semiflexible accuracy (`<= 1.03%`).
  - `rel_freedman2019_crosslinker_stretching_eqS8` (`same_as`): `freedman2019/crosslinker_potential_si_eq4` (`U_stretch`) ↔ MEDYAN Eq. (S8).
- **Parameter contrasts & validations (`q_freedman2019_medyan_vs_afines_comparison_and_no_validation`)**: AFINES softens actin stretching **20,000×** (`5 pN/um` vs `100 pN/nm`) and crosslinker stretching **800×** (`10 pN/um` vs `8 pN/nm`) with **9.26×** longer segments (`1 um` vs `108 nm`), quadratic excluded volume, and force-independent two-step MC unbinding (`0.05 s^-1`). 0 cross-validations (`freedman2019`'s 8 `origin: measured` datums are 2 TIRF actin polymerization speeds and 6 passive fascin/α-actinin sorted domain lengths, none predicted by `chandrasekaran2019`'s runnable models).
