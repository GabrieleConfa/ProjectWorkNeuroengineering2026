# Design of the 3D simulation and the viewer — status and next steps

Status: **2026-11-26** — Gmsh model (D-14), **all 48 contacts solved** for UA01, viewer data exported
for all of them; plan for configuration selection and block agreed (§5.6, D-15). Earlier: 2026-11-22
(oval nerve and rigid tangent cuff pads, D-10/D-11); open questions updated 2026-11-23; roadmap from
the 2026-11-23 tutor meeting in §5. One page to know what exists, how the pieces fit, and what to do next.
The details are elsewhere: model and code walkthrough in [`model_3d.md`](model_3d.md), choices in
[`decisions.md`](decisions.md) (D-NN), numbers and sources in [`parameters.md`](parameters.md), how to
run things in the [README](../README.md).

## 1. The pieces and the data flow

```
config/pw6_3d_params.m          every parameter and modelling choice
        │
        ▼
src/main/main_pipeline_3d.m     patient → 3D geometry → mesh → FE solve (1 A per contact) → figures
        │  writes
        ├── data/patients/<id>.mat            the virtual patient (versioned)
        ├── results/main3d/<id>_3d.mat        mesh + potentials, 80 MB for 6 contacts (local only)
        └── results/main3d/<id>_3d_*.png      figures (versioned)
        │
        ▼
scripts/export_viewer_data.m    samples V, |E|, |J|, AF on a regular grid inside the nerve
        │  writes
        └── results/viewer/<id>/geometry.js + AS<k>.js     ≈ 7 MB per contact (versioned)
        │
        ▼
viewer/index.html               three.js viewer, opens from disk (double-click), offline
```

Why each step is separate:

- **The solve is the expensive part** (≈ 3–3.5 min and several GB per contact) and is done once; the
  problem is linear, so every current is a rescaling of the 1 A solution (D-07).
- **The export** turns the unstructured FE solution into regular grids that a browser can render as
  3D textures. |E| and |J| are taken from the element containing each grid point (exact P1
  gradient), not from finite differences on the grid, which would turn the potential jump across the
  perineurium (D-06) into fake bright rings. Locating the grid in the mesh takes ≈ 4 s
  (`locateRegularGrid`, was ≈ 15 min with `pointLocation`) and is cached per run.
- **The viewer** only displays: no physics runs in the browser. Its data files are JavaScript
  because a page opened from disk may load scripts but not fetch files.

## 2. What is done

| Part | State | Where |
|---|---|---|
| Virtual patient (10 fascicles, mid-upper arm, literature-based, seeded) | done, tested | `src/geometry/`, D-09 |
| 3D geometry (70 mm oval nerve, 50 mm insulating cuff shaped on it, 48 contacts; rigid cuff pads tangent to the oval) | done | `src/mesh/buildNerveGeometry3D.m`, `src/electrode/placeCuffContacts.m`, D-02…D-05, D-10, D-11 |
| FE assembly with perineurium thin layer | done, exact-solution tests | `src/fem/`, `tests/testThinLayerFEM.m`, D-06 |
| Solve, checks (current conservation, equipotential contacts, access resistance) | done for **all 48 contacts** (Gmsh mesh, 2026-11-26; ≈ 30 s per contact with pcg) | `main_pipeline_3d.m`, `solve3DContacts.m`, `model_3d.md` §10 |
| Browser viewer: rotate, stimulate (Space/click), V, \|E\|, \|J\|, AF, transverse and longitudinal cuts | done, data versioned | `viewer/`, `model_3d.md` §8.2 |
| Activating function on the grid + "AF reached by 10% of the fibres" per fascicle | prototype in the export (printout only) | `export_viewer_data.m` |
| Phase 5 proper (fibres, diameters, AF threshold, thresholds in µA/nC) | **not started** (skeletons) | `src/activation/` |
| Safety (Shannon, CIC) | **not started** (skeleton) | `src/safety/` |

Result with the oval nerve (UA01, station −4, AF reached by 10% of the fibres, V/m² per µA): a
contact inside a fascicle drives that fascicle far more than the others (AS 2 in F7: 58×; AS 1 in F2:
23×; AS 4 in F1: 19×); a contact across a fascicle boundary drives two fascicles (AS 3: ratio 2.1);
the cuff pads drive the superficial fascicle under them (2–3.5×). Table in `model_3d.md` §8.2 (the
circle results of 2026-11-20 are kept there for comparison).

## 3. Next session: the full simulation (all 48 contacts) and heatmap B

**Status 2026-11-26: steps 1–3 done.** Every contact is saved as soon as it is solved and reloaded on
the next run (`solve3DContacts.m`); all 48 contacts were solved on the Gmsh mesh with pcg (≈ 30 s each);
the viewer data were re-exported (≈ 2 MB per contact, gzip, ≈ 91 MB in total, committed). Lesson from
that run: the figures of `main_pipeline_3d` §5 took hours because every contact located its sampling
points in the 4.5 M-element mesh again with `pointLocation`; they are now located once with
`locateRegularGrid` (≈ 1 min for the 8 figures), and `<id>_3d.mat` is saved before the figures
(`model_3d.md` §5 walkthrough). Step 4 (heatmap B) is superseded by the MRG threshold table of §5.6 (step A1).

Goal: the **contact × fascicle heatmap of activation strength** ("heatmap B"): for each of the 48
contacts and each fascicle, the activating function reached by 10% of the fascicle's fibres,
normalised (e.g. per contact, or as relative threshold 1/AF). It shows which contact reaches which
fascicle and whether every fascicle is reached by at least one contact — the project question —
without needing the AF threshold value yet (that turns it into thresholds in µA later).

Steps:

1. **Make the solve loop robust before the long run** (`main_pipeline_3d.m` §4). The first run nearly
   ran out of memory (31 GB laptop, ≈ 665k unknowns, direct factorisation per contact); the oval mesh
   is larger (863k nodes, 4.56 M tetrahedra) but the six contacts of station −4 solved fine on
   2026-11-22 (block-wise assembly now bounds the assembly memory). Where the memory goes: the
   stored solutions are small (48 contacts × 863k nodes in single precision ≈ 170 MB); the peak is
   the direct (Cholesky) factorisation of one contact's matrix, whose fill-in takes many GB, plus
   the temporary copies of `K` and `K(free,free)`. Saving per contact is therefore for crash safety
   and resuming, not for RAM. Changes still to make:
   save each contact as soon as it is solved (one file per contact, or append) and skip contacts
   already on disk, so a crash or a stop loses nothing and the run can resume; free the per-contact
   matrices between solves; optionally try `p.solve3d.solver = "pcg"` on one contact and compare with
   the direct solution (time, memory, difference in V).
2. **Run all 48 contacts**: `p.solve3d.sections = []` in `config/pw6_3d_params.m`; ≈ 48 × 3.5 min ≈ 3 h
   plus geometry and mesh, best overnight. Alternative to decide: solve only the 24 contacts of
   stations −4…−1 and use the mirror symmetry z → −z of the layout for the other 24 (half the time,
   but the mesh is not exactly symmetric, so the mirrored fields are an approximation).
   Note: the mesh is not identical between runs (446k vs 510k nodes seen), so the six station −4
   results will change slightly; the whole set must come from one run.
3. **Re-export for the viewer**: `export_viewer_data` (new mesh → the location cache is rebuilt,
   a few seconds since `locateRegularGrid`). 48 contacts × ≈ 7 MB ≈ 320 MB of viewer data: decide before committing whether to
   commit all of it, coarsen the grid, or commit only the final run.
4. **Heatmap B**: move the "AF reached by 10% of the fibres" computation from the export printout
   into a function of `src/activation/` (WP4's module: agree with its owner) that takes the solution
   and returns a 48 × nFascicles table; plot it as a heatmap (contacts grouped by station and type,
   fascicles ordered by size or position); save to `results/`. Open choices for it, to log in
   `decisions.md`: fibres on a regular grid vs random positions (`makeFibres`), one fibre diameter
   (Δz = 1 mm ↔ 10 µm) vs a diameter distribution, normalisation of the colour scale.
5. Then phase 5 proper: AF threshold per fibre diameter (candidate source: Peterson 2011, not yet in
   `Papers/`) → threshold currents and charges, validation against 5.39 ± 0.98 nC (Zelechowski 2020),
   safety check.

## 4. Open questions

- Tutor: shank body size and material, cuff wall thickness, does the shank pierce or push aside the
  fascicles (D-09 assumes pierce), stimulation pulse width, orientation of the device on the oval
  nerve and closing force of the cuff (D-10, D-11) — list in `project_brief.md`. **Partly answered
  at the 2026-11-23 meeting** (probe 100 µm, cuff in PDMS, the probe bends): see §5.
- To verify: roundness 0.62 in Pusnik 2025 Table 1 (taken from Marco's branch; paper not in `Papers/`).
- Convergence study (D-08): Marco's branch has one (`scripts/run_convergence_study.m`, 4 meshes,
  Richardson extrapolation) for his own model; adapting it to `main_pipeline_3d` needs the pipeline
  split into functions — to do before quoting thresholds. Proposed method (2026-11-23, to agree and
  log under D-08):
  - **Levels.** Scale all `p.mesh3d` sizes by one factor s (e.g. 1.6, 1.25, 1, 0.8; ratio ≥ 1.3
    between levels), with 4 levels rather than 3 so that we can see whether the trend is monotonic.
    For an unstructured mesh, use an effective size h = (volume / number of elements)^(1/3).
  - **Remeshing noise.** Our meshes are not nested: every level is a new tetrahedralisation, and
    the same geometry already gave 446k vs 510k nodes. So mesh the reference level 2–3 times and
    take the spread as the noise floor. A change between levels below it says nothing.
  - **Quantities.** Evaluate them at the same points and fibres on every mesh: V at the active
    contact, V at a near and a far fascicle centre, the 10 %-of-fibres activating function (AF) for
    2–3 contact/fascicle pairs, and the rank correlation of heatmap B between levels.
  - **Error estimate: Richardson in the GCI form, with checks** (Grid Convergence Index; Roache,
    Celik et al. 2008, J. Fluids Eng., not yet in `Papers/`, needs a note in `docs/literature/`).
    From the three finest meshes compute the observed order p. If p is plausible (≈ 0.5–3) and
    the changes are monotonic, report the extrapolated value with its GCI as the uncertainty.
    Otherwise do not extrapolate: report the change between the two finest meshes (the criterion
    of Davis 2023: < 0.5 % on the threshold). Why not plain Richardson: it assumes smooth
    convergence in the asymptotic range, which non-nested meshes, the near-singular contact
    edges, the thin perineurium and the AF (a second difference of a piecewise-linear V) may not
    give.
  - **Per-mesh sanity check.** The current through the ground must equal the injected current
    (already printed as `conservErr_pct`).
  - **To decide.** The tolerance (proposal: 5 % on the AF, 2 % on V). The finest level will not fit
    in 32 GB with the direct solver (memory grows faster than the number of nodes). The proposed way
    round is to refine only the fine zone |z| < 12 mm and/or use `pcg` for that level.
- **Hardware (2026-11-23).** A desktop with 32 GB RAM and a GPU has the same RAM as the laptop, so the
  same limits apply. Its advantages are long unattended runs and probably more CPU cores (the direct
  solver is multithreaded). The GPU is not usable as things stand:
  - the Parallel Computing Toolbox, which MATLAB needs for GPU computing, is not installed on the
    laptop (`gpuDeviceCount` is undefined; check `ver` on the desktop), and it needs an NVIDIA card;
  - sparse direct solves have little GPU support, the `ichol` preconditioner is built on the CPU,
    our matrix is badly conditioned (thin perineurium, contacts at 1e3 S/m) and consumer GPUs are
    slow in double precision;
  - geometry, meshing and grid location are CPU-only.

  Run all 48 contacts on a single machine: the meshes differ between runs.
- **Multi-contact configurations (simultaneous, any currents, monopolar/bipolar/multipolar).**
  - **Superposition.** The problem is linear, so any configuration is post-processing of the 48
    single-contact 1 A solutions: V = Σ Iₖ·Vₖ, and likewise AF = Σ Iₖ·AFₖ. This includes bipolar
    configurations, V = I·(V_a − V_b), which need no new solve: the two returns through the
    grounded boundary cancel (net ground current Σ Iₖ = 0). The grounded outer surface barely
    matters for a balanced pair (the dipole field ∝ 1/r² is ≈ 0 there; see D-04 above).
  - **Voltage-controlled stimulation.** Use the 48 × 48 transfer-resistance matrix R (Rᵢⱼ =
    potential of contact i for 1 A on contact j): solve R·I = V_imposed, then sum the fields.
  - **Limits.** Activation is not linear: sum the fields first, then threshold; never add the
    recruitment of single contacts. Only overlapping pulses add up; interleaved pulses are separate
    stimuli.
  - **Open: idle contacts "insulating" (D-07, as in 2D) vs "floating"** (`p.electrode.idleContacts`).
    Superposition is exact only if the system matrix is the same for every contact. With insulating
    idle contacts it is not: in Vₐ contact b is an insulator, in an a+b configuration it is a
    conductor. The sum is then an approximation, worst for neighbouring contacts on the same shank,
    which are the typical bipolar pairs.
    - With "floating", every contact is always a conductor, so superposition is exact and the matrix
      can be factorised **once** for all 48 contacts (≈ 10–15 min instead of ≈ 3 h; needs a code
      change, today it re-factorises even with "floating"). It is also arguably more physical: a
      disconnected metal contact is equipotential with zero net current.
    - Cost of "floating": it departs from the 2D model and changes D-07, and it ignores the
      electrode–tissue interface impedance, so it slightly overestimates the shunting by idle
      contacts.
    - **Proposed test:** solve one contact (AS 2) both ways and compare V and the 10 %-of-fibres AF
      per fascicle. Decision for the group (and possibly the tutor).
    - If we keep "insulating", check the configurations that matter with one dedicated solve each
      (≈ 3.5 min).
  - A new solve is needed only if the matrix changes: contact states (above), tissues or geometry
    (new patient, new mesh), or non-linear interface physics (not in the model).
  - **Free verification test** once all 48 contacts are solved: R must be symmetric (Rᵢⱼ = Rⱼᵢ,
    reciprocity). Add it to `tests/`.
- Group: 24 contacts + symmetry or 48 contacts; how much viewer data to version; who owns heatmap B
  (WP4 activation module).
- Known numerical points not yet checked: mesh convergence (D-08), ground vs insulated ends (D-04),
  location of the "poor shape quality" elements reported by the mesher.

## 5. Roadmap after the tutor meeting (2026-11-23)

What the meeting asked for is written as stated; "Notes" are our own remarks (to check, not tutor
statements). Order: items 1–3 are for the next session, 4–7 come later.

### 1. Electrode probe that bends inside the nerve (in the 3D model since 2026-11-26: Gmsh mesher, D-14)
**Status 2026-11-26 (D-13):** Marco's insertion Monte Carlo is on main (`src/insertion/`, study
`scripts/run_insertion_3d.m` on our patient: his force model, and the tutor's rule as a probability
P(radius, angle) = 1/(1+(r50/r)^n)·cos^m(γ), `src/insertion/pierceProbability.m`). It is
2D per station. **D-14:** the default 3D mesher is now Marco's Gmsh model, which builds the bent PDMS lead,
the displaced and irregular fascicles and the PDMS sleeve from that insertion (one draw per station),
with our solver. Still to do: the viewer (circles, no lead body), statistics over several insertion
seeds and patients.

Model as described at the meeting:
- the probe (shank) enters the nerve and, at each fascicle on its way, pierces it with a
  probability that grows with the fascicle radius (a large fascicle is more likely pierced, a small
  one pushed aside);
- the exit point is on the same straight line as the entry (entry and exit stay where the straight
  probe would be);
- between them the probe bends according to the position of the fascicles it does not pierce;
- the probe is **100 µm**; the cuff is **PDMS, high impedance**.

Notes:
- **This reverses current decisions:**
  - D-09: "pierce" always;
  - D-11: no shank deflection. The group had rejected the "deflect" scenario on 2026-11-22.

  Both must be revised in `decisions.md` when this is implemented.
- **Start from Marco's branch** `cuff-oval-viewer-convergence`, which already has a shank "deflect"
  scenario: `docs/insertion_evidence.md`, `config/quicklook3DParams.m`, `scripts/run_3d_quicklook.m`.
- **The PDMS cuff is consistent with the model:** the wall is already a perfect insulator (void,
  D-05).
- **To ask the tutor:**
  - is 100 µm the width or the thickness of the probe?
  - the probe material;
  - the law for the pierce probability: linear in r? relative to which radius? from a paper?
- **Consequences for the model:**
  - the shank body must enter the model as an insulator: today there are contacts only (D-07);
  - the contacts follow the bent probe (position and orientation);
  - the bend becomes part of the virtual patient (random, seeded: reproducible from
    `data/patients/`);
  - to decide: does the bend change from station to station along z?

### 2. Mesh-convergence algorithm (written 2026-11-26, to run)
Method proposed in §4 (GCI with checks) + Marco's algorithm: `scripts/run_convergence_3d.m`, on the
pipeline split into `src/main/build3DModel.m` and `solve3DContacts.m` (D-13). Also measures the
remeshing noise floor and the effect of 5 µm edge refinement (`p.mesh3d.hEdge_m`, off in production).

### 3. Blocking the signal towards the amputated limb (user's part, next session)
**Status 2026-11-25: first round implemented** on the current model (no new solve): the two sites
are the two mirror halves of the layout, distal contacts by z-mirror, MRG fibres in NEURON.
Walkthrough and results: `docs/block.md`; choices: D-12. The "all-round probes" below turned out to
be the existing contacts (group, 2026-11-25). **Next steps (2026-11-26): task B of the plan in §5.6**
(re-run on the Gmsh solution, then one-directional distal stimulation with the configuration optimiser).

Clinical scenario (amputee, after targeted muscle reinnervation — the surgeon reconnects the nerve
to other muscles):
- the electrode has two macro sites, i.e. two all-round (circumferential) electrode probes;
- the section closer to the shoulder (proximal) is stimulated to give the sensation;
- the distal section, near the severed limb, is driven with a reversed anode–cathode configuration
  to stop the action potentials travelling towards the reconnected muscles, so the patient feels no
  discomfort from them.

Notes:
- **Where this sits in the literature.** It is unidirectional stimulation / conduction block:
  - anodal block with an asymmetric configuration (classic reference: Sweeney & Mortimer 1986,
    not in `Papers/`);
  - alternatives: kilohertz high-frequency (KHFAC) block, collision block.
- **The static FEM is not enough.** It gives V and the activating function: positive =
  depolarisation, negative = hyperpolarisation, the region where a block could occur. Whether an
  action potential is actually blocked needs a dynamic fibre model (e.g. MRG). Marco's branch has
  `scripts/run_mrg_fibre.py`.
- **Anodal block needs high currents.** That brings the safety check (Shannon, CIC) and the risk of
  exciting at the edge of the hyperpolarised region (virtual cathode, anode break).
- **To ask the tutor:** geometry of the two macro sites / all-round probes (size, position along z,
  material).

### 4. Non-linear electrode–tissue interface (later)
Notes:
- **It is a contact boundary condition, not a tissue property.** The interface is the double layer
  (capacitance + faradaic current, e.g. Butler–Volmer).
- **It breaks linearity.** Superposition (§4) no longer holds, so every configuration and amplitude
  needs its own, time-domain solve.

### 5. Frequency-dependent (and field-dependent) tissue model — plan of 2026-11-26
Status: **plan, not implemented, not yet decided** (it gets a decision number when adopted). The user proposed
div(σ(|E|, ω) grad V) = 0. Plan: the ω part as a linear capacitive model (below); the |E| part
not in the FEM but as a post-processing safety map, already written (`scripts/run_field_safety.m`).

**Formulation.** ω stays real; what becomes complex is the admittivity σ*(ω) = σ + jωε₀ε_r of each
tissue (anisotropic σ as now), and the potential becomes a phasor V̂:
div(σ*(ω) grad V̂) = 0. Still quasi-static (no induction) and still **linear**, so superposition
holds. Equivalent time-domain form: div(σ grad V) + div(ε grad ∂V/∂t) = 0. The perineurium thin
layer becomes a parallel RC sheet: J = G_s ΔV + C_s ∂ΔV/∂t, with C_s = ε₀ε_r / t.

**Why it matters (sources in `docs/literature/`).**
- Perineurium (Horn 2023): ε_r ≈ 2018 → C_s = 8.85e-12 × 2018 / 13.8e-6 ≈ 1.3e-3 F/m²
  (0.13 µF/cm², our recomputation; the paper's 13.5 µF/cm² is the open issue H1). With **our**
  sheet resistance (1149 Ω·m × 14 µm = 161 Ω·cm²) τ = R_s·C_s ≈ 21 µs (corner ≈ 7.7 kHz); with
  Horn's 510 Ω·cm², τ ≈ 65 µs (2.4 kHz). Either way τ is of the order of the 50 µs stimulation pulse.
- Epineurium (Horn 2023): ε_r ≈ 9.4e6 (± 8.2e6). With our σ_epi = 0.159 S/m the corner is
  σ/(2πε₀ε_r) ≈ 300 Hz: above it the epineurium is dominated by its capacitance, i.e. during a
  short pulse it conducts much better than in the resistive model. A constant ε_r overstates this at
  high frequency; Lazorchak 2022 used low- and high-frequency resistivities (epineurium ρ₀ = 68,
  ρ∞ = 1.3 Ω·m; perineurium ρ₀ = 3460, ρ∞ = 273 Ω·m) with constant ε_r, and the best fit to their
  in-vivo thresholds was **low-frequency ρ + measured ε_r** (parallel RC per tissue). We follow that.
- Lazorchak 2022 (time-dependent FEM + MRG): permittivity matters mainly at short pulse widths.
- **Kilohertz (KHFAC) block**, the standard alternative to anodal block, needs this model.

**Missing parameters (to find before implementing).** Endoneurium ε_r (Horn does not measure it);
muscle ε_r (candidate: Gabriel, Lau & Gabriel 1996, Phys Med Biol 41:2271, parametric models, also
in the IT'IS tissue database; not in `Papers/`). Do the course slides give any permittivity?

**Implementation (two forms, both on our own P1 assembly).**
1. **Transient, for pulses.** Implicit Euler: (K_σ + K_ε/Δt) Vⁿ⁺¹ = (K_ε/Δt) Vⁿ + Fⁿ⁺¹, where K_ε is
   assembled exactly like K_σ (`assembleConductionP1` with ε instead of σ) and the thin layer adds
   G_s + C_s/Δt. The matrix is **real, symmetric positive definite and constant**: factorise once per
   contact (`decomposition`), then each step is a pair of triangular solves; memory as now. Output:
   the **step response** of each contact at the fibre samples. Two step sizes (e.g. 5 µs up to
   ~0.2 ms, then 50 µs up to ~3 ms, since the epineurium τ = ε/σ ≈ 0.5 ms) → ~100 steps and two
   factorisations per contact. Option `p.solve3d.mode = "static" | "transient"` in
   `main_pipeline_3d`, cached per contact like the static solve.
2. **Phasor, for KHFAC.** One complex solve per contact at the chosen frequency (complex symmetric,
   ~2× the memory of the real factorisation: check RAM first).

**Coupling to the fibre model.** With capacitance, the potential of a current pulse changes shape in
time, so I(t)·V(z) no longer holds. But the system is linear and time-invariant, and our pulses are
sums of steps (rectangles, staircases), so V(z, t) = Σ_c Σ_k ΔI_c,k · s_c(z, t − t_k), with s_c the
step response of contact c. `mrg_fibre` then applies a time-varying `e_extracellular`
(`Vector.play` per compartment). Storage: ~30 log-spaced times × 453 fibres × 1321 samples ×
4 bytes ≈ 70 MB per contact.

**Step 0: decide whether it matters.** Transient solve of AS 2 only; recompute the F7 thresholds
(50 µs) and the first-round anode runs (onset transients, virtual cathodes). If the 10 % threshold
moves by more than ~10 %, adopt the transient model for everything; otherwise document that the
quasi-static assumption is justified for our pulses (a result for the report) and keep the phasor
form only for KHFAC.

**Field dependence σ(|E|): not in the FEM, a safety map instead.** The physical mechanism is
electroporation (membrane pores at fields of the order of hundreds of V/cm; threshold for nerve
tissue still to be sourced: `p.safety.electroporation_Vpm = NaN`). Rough estimate at the contact
surface (≈ 0.12 mm², σ_endo ≈ 0.17 S/m): ≈ 5 V/cm per 10 µA, i.e. negligible at stimulation
currents but ≈ 500–2500 V/cm at the 1–5 mA of the block sweeps, within ~100 µm of the contacts. So
the nonlinearity appears only where the tissue is damaged anyway, while modelling it would remove
superposition (one 3D Newton solve per pattern and current) and would still be crude (pores open and
reseal in time). `scripts/run_field_safety.m` (functions `src/safety/elementFieldP1.m`,
`fieldExposure.m`, `exposedVolume.m`, test `tests/testFieldSafety.m`) computes, for each pattern,
the volume of tissue above several field levels vs current and the current at which 1e-3 mm³
exceeds each level (the peak field at the contact edges is mesh-dependent, the volume is not).

**Two different frequencies.** "Frequency" here is the spectral content of one pulse (kHz) or the
KHFAC carrier, not the pulse repetition rate of item 6.

### 6. Choosing the configuration to reach one fascicle (the main goal)
Given a target fascicle (connected to a given muscle or skin territory), decide which contacts to
turn on, with which currents, pulse width and rate. Example: sensory feedback from a prosthesis that
feels like the hand closing, which needs specific fascicles.

Notes:
- **An optimisation problem (inverse problem / input allocation).** The activating function is
  linear in the contact currents (§4), so: maximise recruitment of the target, keep the other
  fascicles below threshold, respect the safety limits per contact (charge, charge density).
  These are linear constraints, i.e. a linear/quadratic program.
- **Needs:**
  - the AF thresholds per fibre diameter (phase 5);
  - a map fascicle → function (which muscle or skin area each fascicle serves). In the upper arm
    the fascicles exchange fibres along the nerve, so that map needs a source (to verify).
- **What each parameter controls.** Amplitude and pulse width set how many fibres are recruited;
  the pulse rate sets the firing rate, i.e. the perceived intensity (rate coding). The FEM solution
  per pulse does not depend on the rate; the rate matters through fibre dynamics.

#### Plan agreed on 2026-11-26 (configuration selection + block; choice record D-15)

**One idea for both tasks.** The model is linear, so for any configuration V = Σ Iₖ·Vₖ over the 48
unit solutions, and along every fibre the activating function (AF) is the same weighted sum of 48
per-contact profiles. Searching over configurations is therefore cheap with the AF; every candidate is
then **verified with the MRG fibres in NEURON** (`src/block/mrg_fibre.py`), because the AF is only a
proxy (the true threshold also depends on the field shape and the pulse width).

**The device works at fascicle level (group, 2026-11-26).** It cannot stimulate or record single
fibres: it sets contact currents and can only observe whole fascicles (compound activity). Consequences:
- the fibres in the model are a **means**, not an output: each fascicle is a population of model
  fibres, and every result is reported per fascicle (fraction of its fibres recruited / blocked,
  10%-of-fibres criterion as in Zelechowski 2020);
- the optimiser's targets and constraints are **fascicle recruitment fractions** (e.g. target ≥ 50%,
  every other fascicle ≤ 10%), not individual fibres;
- physically the stimulus still acts fibre by fibre (the fibres closest to the contact fire first), so
  a partial recruitment inside a fascicle can happen, but the device cannot measure or exploit it:
  a function whose fibres are mixed with others inside a fascicle (Delgado-Martinez 2016,
  `src/patient/functionCatalog.m`) cannot be separated from them; the best possible is to choose the
  fascicles where that function's share is largest;
- recording (§7) is a fascicle-level signal too: by reciprocity the 48 unit solutions are the lead
  fields of the compound activity of each fascicle.

**Task A — configuration that fires a chosen fascicle**
- **A1. Monopolar screening with MRG.** For every contact and every fibre, the threshold current by
  bisection → recruitment curve of every fascicle → contact × fascicle table of thresholds (µA, nC),
  selectivity index and selectivity window (how far the current can rise before the next fascicle
  reaches 10%; AS 2 → F7 gave ×12 in the first block round). Cost ≈ 48 × 453 fibres × ~10 runs × 0.3 s
  / 7 processes ≈ 2–3 h per fibre diameter (overnight). It also **calibrates our AF threshold** (AF at
  the MRG threshold), which removes the need for an external AF-threshold source (Peterson 2011). It
  supersedes heatmap B (§3).
- **A2. Multi-contact current steering as a linear program** (`scipy.optimize.linprog`, Python, next to
  NEURON). Variables: signed currents Iₖ. Linear in Iₖ: AF ≥ θ on the target fascicle's fibres
  (maximise the margin or the recruited fraction), AF ≤ θ/(safety factor) on the other fascicles' fibres
  at every node, |Iₖ|·PW ≤ the charge limit per contact (CIC of the contact material, Shannon), at most
  N active contacts, optionally Σ Iₖ = 0 (local return). Then MRG verification, adjust θ, 2–3 loops.
- **A3. Target a function** (group of muscles or skin area) instead of a fascicle: target set = the
  fascicles with the largest share of that group (`src/patient/`), within the fascicle-level limit above.
- **A4. Robustness**: A1–A2 over several lead insertions and virtual patients (phase 6): how often each
  target can be reached selectively, not only the best case on UA01.
- **A5. Viewer**: a "target" drop-down that loads the optimised pattern into Pattern mode, through the
  configuration JSON shared with NEURON (§8).

**Task B — block towards the amputated limb** (first round and its results: `block.md` §6, §5.3 here)
- **Update 2026-11-27 (tutor, D-16):** the scheme to test first is the **anodal block on the same
  distal contacts**, judged per fascicle (residual volley at the muscle < 10 % of the target fascicle,
  no new activity in the others). Written: `notebooks/block_volley.ipynb` (includes B0 check and B2
  ramps). B1 stays as the alternative if the anode cannot meet the criteria.
- **B0.** Re-run the first round on the Gmsh solution (all 48 contacts solved, no mirror), ≈ 10 min.
- **B1. One-directional distal stimulation (collision block, the most promising).** A distal cathode
  with anodes on its distal side (e.g. cathode on station +3, anodes on +4): the distal AP travels only
  towards the brain and annihilates with the descending AP, so nothing reaches the muscle and the brain
  receives only the proximal AP. It is task A with extra constraints: the distal cathode must fire the
  same fascicle the proximal site fired and nothing else (otherwise extra sensations), the distal anode
  must hyperpolarise it enough to stop the downward branch without virtual-cathode firing (anode spread
  over several contacts or the cuff pads = lower charge density). Timing window = travel time between
  the sites (≈ 0.4 ms over 20 mm at 52 m/s for 10 µm fibres; wider for smaller, slower fibres).
- **B2. Pulse shaping**: slowly rising anodal phase (the virtual cathodes accommodate instead of
  firing), slowly decaying end (quasi-trapezoidal, no anode break). Needs piecewise-linear waveforms
  in `mrg_fibre.Channel` (today rectangular only).
- **B3. Fibre-diameter sweep**: excitation and block scale differently with D.
- **B4. Kilohertz block**: later, only as a discussed alternative (needs §5.5; onset response).
- **Safety on every candidate**: charge per phase and charge density vs the contact material and the
  Shannon criterion — in the first round this ruled out most single-contact block settings.

**Order**: (1) B0 → (2) export fibre profiles for all 48 contacts (`export_block_profiles`) → (3) A1
(overnight) → (4) A2 + MRG loop (Python module + notebook) → (5) B1 reusing the A2 optimiser → (6) B2/B3
→ (7) A3/A4 (phase 6) → (8) viewer drop-down + shared JSON.

**Still open before A1** (user / tutor; to log in D-15 when decided):
1. fibre diameters: one 10 µm (as so far) or a distribution per fascicle (needs a median-nerve source);
2. selectivity metric: 10%-of-fibres criterion (Zelechowski 2020) and which selectivity-index
   definition (source to find before citing);
3. optimiser limits: maximum number of active contacts, return to ground or balanced local return,
   contact material for the charge limit;
4. block goal (tutor): block only the fibres the proximal site fired, or the whole nerve; **answered
   2026-11-27 (D-16)**: anode on the same distal contacts, the global activation of the target
   fascicle must stop (partial block allowed); notebook `block_volley.ipynb` (`block.md` §5);
5. **what the stimulator can drive** (user, 2026-11-26: it may control the electrode only as a whole,
   not contact by contact): single contacts or only groups, opposite polarities on nearby contacts,
   number of independent channels, return electrode, pulse shapes. Questions in `project_brief.md`
   ("How the stimulator drives the contacts"). A1 (monopolar screening) is useful in any case; A2, B1
   and B2 must be re-scoped to what the hardware allows once this is answered.

### 7. Recording and stimulation are time-multiplexed
The device cannot record and stimulate at the same time (the stimulus artefact saturates the
recording amplifiers), so the stimulation pulses are scheduled in the gaps between recording
windows.

Notes (to verify):
- **The two rates are different things.** The stimulation pulse rate (typically tens to hundreds of
  Hz) is not the sampling rate of the recording.
- **The sampling rate.** Neural recording is usually sampled at about 20–30 kHz per channel; the
  "up to 100 kHz" from the meeting may be the device total, or a specific device.
- **They are linked only by scheduling.** The pulse (width ≈ 0.05–0.3 ms) plus the blanking time
  must fit in the gaps.
- **For recording, the model is already usable.** By reciprocity, the sensitivity of contact i to a
  fibre source at point x equals the potential at x when 1 A is injected into contact i. So the 48
  stored solutions are also the "lead fields" for simulating recording.

### 8. End goal: real-time simulation in the viewer, any electrode configuration (user, 2026-11-25)
**Status 2026-11-26:** first part done — the viewer has a *Pattern* mode: any set of contacts with
signed currents, summed in the browser from a 16-bit potential exported per contact (gzip), with V,
|E| (finite differences), AF+ and AF− (README, *Patterns*). **2026-11-27:** first step of the fibre part —
*Volley* playback of NEURON runs (AP times at every node, exported by `volley.export_viewer_volley`;
README, *Volley playback*); the exported runs will also be the reference for the JS port. Still to do for this goal: the fibre model in
the browser (JS port of MRG checked against NEURON), time courses (pulse widths, delays), and the
"floating" idle contacts for an exact superposition (today: idle contacts insulating, D-07).

The final viewer must let the user set any configuration (which contacts, currents, polarities,
pulse timing on each site) and see the result in real time. What this implies:

- **Fields: real time by superposition.** V for any configuration is Σ Iₖ·Vₖ of the 48 unit
  solutions, so the browser only sums precomputed volumes (a GPU shader) and derives E and the
  activating function from the summed V. Consequences:
  - all 48 contacts must be solved (or 24 + z-mirror, as the block study already does, D-12);
  - the per-contact V must be shipped **linearly coded with enough precision to be summed**
    (16-bit, not today's 8-bit), and E, J, AF computed in the viewer, not precomputed per contact
    (the log-coded E/J cannot be summed); size ≈ 1.6 M voxels × 2 B per contact (≈ 80 MB for 24 +
    mirror) — to check against the WebGL texture limits, or coarsen the grid;
  - superposition must be exact: this is a strong argument for **"floating" idle contacts** (§4,
    multi-contact configurations), which also allows one factorisation for all 48 contacts.
- **Fibre responses: the NEURON model cannot run in the browser.** Plan: Python/NEURON (MRG,
  `src/block/mrg_fibre.py`) stays the **reference**; the viewer gets a JavaScript port of the same
  fibre model (MRG double cable, or a node-only reduction if MRG is too slow), checked against the
  NEURON results (thresholds, AP times) in a test. Rough cost of MRG: ~570 compartments × 800 time
  steps per fibre → milliseconds per fibre in a Web Worker, so ~100–450 fibres in about a second.
  The viewer would then animate the APs along the fibres (the raster of `docs/block.md` in 3D).
- **One configuration format for both sides.** A configuration = a list of channels, each = contact
  weights (sign = polarity) + rectangular pulses (t_on, t_off, amplitude). `mrg_fibre.Channel` already
  uses this; the viewer should read and write the same JSON, so that a configuration found in the
  viewer can be re-run in NEURON and vice versa.
- **Fibre data for the viewer.** The per-contact profiles along the fibres
  (`results/block/<id>_fibre_profiles.h5`: 453 fibres × 1321 samples × contacts) are too large to
  ship as they are; a coarser set (fewer fibres, 100 µm samples, 16-bit) will be exported for the
  viewer.
