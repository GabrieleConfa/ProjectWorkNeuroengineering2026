# Work plan and division of work

**Status: proposal (2026-11-08), to be agreed at the kick-off meeting (2026-11-10).**
Group of 5. Interfaces between packages: `docs/architecture.md`.

## 1. How the work is split

Five **work packages (WP)**, one owner each. The packages follow the modules of the code, so each
person owns a folder of `src/` and nobody waits for someone else's file to change. To keep the
load even over time, every package has an **early part** (needed for the first 3D simulation) and
a **late part** (needed for the results), and every owner also writes the **report section and
the slides** of their package.

| WP | Title | Owner | Folders | Early part (sprints 1–2) | Late part (sprints 3–4) |
|---|---|---|---|---|---|
| WP1 | Anatomy and virtual patients | _(name)_ | `src/geometry`, `data/patients` | tutor's section → `makePatient`/`buildCrossSection2D` (patient P00); generator design from Delgado-Martínez 2016 | random cross-section generator (seeded, no overlaps, validated against the paper's statistics); population of virtual patients |
| WP2 | Electrode, 3D geometry and mesh | _(name)_ | `src/electrode`, `src/mesh` | extrusion to 3D, insertion of the 48 AS (+ shank?), 3D plot; mesh sizing | mesh-convergence study; perineurium representation in the mesh; run-time budget per patient |
| WP3 | FEM solver and verification | _(name)_ | `src/fem`, `tests` | `createpde` set-up with anisotropic σ, boundary conditions, current injection on the AS, unit-current solve | verification suite (point source, linearity, symmetry ±z, convergence); batch solve of all AS × patients |
| WP4 | Fibres, activation and thresholds | _(name)_ | `src/activation` | fibre model; activating function tested on the analytical point source; literature for the AF threshold and fibre diameters | thresholds and recruitment per fascicle; validation vs 5.39 ± 0.98 nC (Zelechowski 2020); which fascicles each AS reaches |
| WP5 | Tissue properties, safety, integration | _(name)_ | `src/materials`, `src/safety`, `config`, `scripts` | conductivity set and perineurium ρ·t (literature, proposal); `config/default_params.m`; safety literature (Shannon k, CIC of gold) | safety check; sensitivity analysis on conductivities; `run_pipeline`, result figures; report assembly |

Cross-cutting duties, also split evenly:

- **Code review in a ring:** WP1 reviews WP2, WP2 → WP3, WP3 → WP4, WP4 → WP5, WP5 → WP1.
  Nobody merges their own code unreviewed. This also means every module is known by two people.
- **Documentation:** each owner keeps their rows of `parameters.md` and their entries in
  `decisions.md` current. WP5 checks consistency before each tutor meeting.
- **Tutor contact:** one person collects the questions (`docs/project_brief.md` → "Still
  unknown") and sends them; rotate every meeting.

## 2. Sprints

Dates of the tutor meetings and the final deadline are not known yet (agenda item). Sprints are
therefore listed by goal; fix the dates on Tuesday.

### Sprint 1 — "3D nerve with the electrode and a first solution" (goal for the next tutor meeting)

The tutor asked to add the electrode's active sites to the model and to simulate them.

| WP | Task | Done when |
|---|---|---|
| WP1 | Move the tutor's 2D section (7 fixed fascicles, R = 1.5 mm) into `makePatient` (patient P00) + `buildCrossSection2D` | P00 saved in `data/patients/`, face → tissue map tested (the tutor's check of section 9 as a unit test) |
| WP2 | `buildNerveGeometry3D` (extrude P00 + outer domain) and `addActiveSites` (48 boxes); 3D figure of nerve + AS | figure of the 3D nerve with all 48 AS; cuff-site insertion tested |
| WP3 | `solveUnitCurrent` on a toy geometry, then on WP2's; first unit test: point source vs V = I/(4πσr) | V map for one AS on P00 (cross-section at the AS and profile along z) |
| WP4 | `makeFibres` + `activatingFunction` working on the analytical point-source potential; search the AF-threshold source | plot of AF along a fibre vs distance from the source; threshold source proposed |
| WP5 | Conductivity proposal (Davis 2023 used by the tutor vs Zelechowski 2020 vs Horn 2023), perineurium ρ·t table; `config/default_params.m` filled; gold CIC and Shannon k found | `parameters.md` complete for phase 2 and safety; config loads and is used by the others |

Meeting deliverable: the 3D figure (WP2) + the potential of one AS (WP3) + the parameter table
(WP5) + the plan of WP1/WP4.

### Sprint 2 — "Full pipeline on one patient"
All 48 AS solved on P00; mesh convergence; first thresholds per fascicle; first safety numbers.

### Sprint 3 — "Virtual patients and robustness"
Random generator; N patients; batch runs; sensitivity to conductivities and perineurium; which
fascicles stay unreached by any AS.

### Sprint 4 — "Discussion, report, presentation"
Results across patients, limitations (from `decisions.md`), report and slides.

## 3. Collaboration rules (proposal)

- **Shared repository** on GitHub (private, shared at the kick-off meeting), one branch per task,
  merge through a pull request reviewed by the ring partner. Nobody pushes directly to the main
  branch.
- **Same software for everyone:** MATLAB R2025a + PDE Toolbox (the 3D functions used —
  `extrude`, `addCell`, `femodel` — need a recent release; R2023a is the minimum for `femodel`).
- **Large files stay out of git:** `results/` and `Papers/` are ignored; everything in `results/`
  must be regenerable from `scripts/`.
- **Before pushing:** `setup_paths; runtests('tests')` passes.
- **A choice that is not in the brief** → entry in `docs/decisions.md` (proposed by the owner,
  agreed in the group chat or at the weekly meeting).
