# Design of the 3D simulation and the viewer — status and next steps

Status: **2026-11-20**. One page to know what exists, how the pieces fit, and what to do next.
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
        ├── results/main3d/<id>_3d.mat        mesh + potentials, 50 MB for 6 contacts (local only)
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

- **The solve is the expensive part** (≈ 2.5 min and several GB per contact) and is done once; the
  problem is linear, so every current is a rescaling of the 1 A solution (D-07).
- **The export** turns the unstructured FE solution into regular grids that a browser can render as
  3D textures. |E| and |J| are taken from the element containing each grid point (exact P1
  gradient), not from finite differences on the grid, which would turn the potential jump across the
  perineurium (D-06) into fake bright rings. Locating the grid in the mesh takes ≈ 15 min and is
  cached per run.
- **The viewer** only displays: no physics runs in the browser. Its data files are JavaScript
  because a page opened from disk may load scripts but not fetch files.

## 2. What is done

| Part | State | Where |
|---|---|---|
| Virtual patient (10 fascicles, mid-upper arm, literature-based, seeded) | done, tested | `src/geometry/`, D-09 |
| 3D geometry (70 mm nerve, 50 mm insulating cuff, 48 contacts) | done | `src/mesh/buildNerveGeometry3D.m`, D-02…D-05 |
| FE assembly with perineurium thin layer | done, exact-solution tests | `src/fem/`, `tests/testThinLayerFEM.m`, D-06 |
| Solve, checks (current conservation, equipotential contacts, access resistance) | done for **station −4 only (6 of 48 contacts)** | `main_pipeline_3d.m`, `model_3d.md` §10 |
| Browser viewer: rotate, stimulate (Space/click), V, \|E\|, \|J\|, AF, transverse and longitudinal cuts | done, data versioned | `viewer/`, `model_3d.md` §8.2 |
| Activating function on the grid + "AF reached by 10% of the fibres" per fascicle | prototype in the export (printout only) | `export_viewer_data.m` |
| Phase 5 proper (fibres, diameters, AF threshold, thresholds in µA/nC) | **not started** (skeletons) | `src/activation/` |
| Safety (Shannon, CIC) | **not started** (skeleton) | `src/safety/` |

First result (UA01, station −4, AF reached by 10% of the fibres, V/m² per µA): a contact inside a
fascicle drives that fascicle far more than the others (AS 2 in F2: 26×; AS 4 in F1: 15×); a contact
across a fascicle boundary drives two fascicles (AS 3: ratio 1.2); the cuff contacts drive the
superficial fascicle under them (≈ 4×). Table in `model_3d.md` §8.2.

## 3. Next session: the full simulation (all 48 contacts) and heatmap B

Goal: the **contact × fascicle heatmap of activation strength** ("heatmap B"): for each of the 48
contacts and each fascicle, the activating function reached by 10% of the fascicle's fibres,
normalised (e.g. per contact, or as relative threshold 1/AF). It shows which contact reaches which
fascicle and whether every fascicle is reached by at least one contact — the project question —
without needing the AF threshold value yet (that turns it into thresholds in µA later).

Steps:

1. **Make the solve loop robust before the long run** (`main_pipeline_3d.m` §4). The first run nearly
   ran out of memory (31 GB laptop, ≈ 665k unknowns, direct factorisation per contact). Changes:
   save each contact as soon as it is solved (one file per contact, or append) and skip contacts
   already on disk, so a crash or a stop loses nothing and the run can resume; free the per-contact
   matrices between solves; optionally try `p.solve3d.solver = "pcg"` on one contact and compare with
   the direct solution (time, memory, difference in V).
2. **Run all 48 contacts**: `p.solve3d.sections = []` in `config/pw6_3d_params.m`; ≈ 48 × 2.5 min ≈ 2 h
   plus geometry and mesh, best overnight. Alternative to decide: solve only the 24 contacts of
   stations −4…−1 and use the mirror symmetry z → −z of the layout for the other 24 (half the time,
   but the mesh is not exactly symmetric, so the mirrored fields are an approximation).
   Note: the mesh is not identical between runs (446k vs 510k nodes seen), so the six station −4
   results will change slightly; the whole set must come from one run.
3. **Re-export for the viewer**: `export_viewer_data` (new mesh → the location cache is rebuilt,
   ≈ 15 min). 48 contacts × ≈ 7 MB ≈ 320 MB of viewer data: decide before committing whether to
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
  fascicles (D-09 assumes pierce), stimulation pulse width — list in `project_brief.md`.
- Group: 24 contacts + symmetry or 48 contacts; how much viewer data to version; who owns heatmap B
  (WP4 activation module).
- Known numerical points not yet checked: mesh convergence (D-08), ground vs insulated ends (D-04),
  location of the "poor shape quality" elements reported by the mesher.
