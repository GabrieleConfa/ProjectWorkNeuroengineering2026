# Intraneural stimulation of the human median nerve — FEM project

Neuroengineering course project (Politecnico di Milano, A.Y. 2026/2027) — group work,
tutor Luca Sbarsi (luca.sbarsi@polimi.it).

Does the active-site layout of a new intraneural electrode reach **every fascicle** of the human
median nerve? Answered with a 3D finite element model (MATLAB PDE Toolbox), an
activating-function criterion for fibre activation, a safety check, and several virtual patients.

Brief: [docs/project_brief.md](docs/project_brief.md) ·
Architecture: [docs/architecture.md](docs/architecture.md) ·
Work plan: [docs/work_plan.md](docs/work_plan.md) ·
Parameters: [docs/parameters.md](docs/parameters.md) ·
Decisions: [docs/decisions.md](docs/decisions.md) ·
Meetings: [docs/meetings/](docs/meetings/) ·
Paper notes: [docs/literature/](docs/literature/)

**New to the project?** Read `docs/architecture.md` §1–4, find your work package in
`docs/work_plan.md`, open the skeleton functions of your folder in `src/` (their header is the
contract), and run `setup_paths; runtests('tests')`.

## Structure

| Folder | Content |
|---|---|
| `docs/` | brief, architecture, work plan, parameter table, decision log, one note per paper (papers themselves are not in the repo: download them from the DOIs in the notes) |
| `src/geometry` | phase 1 — nerve cross-section generator |
| `src/materials` | phase 2 — conductivity assignment |
| `src/electrode`, `src/mesh` | phase 3 — active sites, 3D geometry and mesh |
| `src/fem` | phase 4 — PDE setup, solution, post-processing |
| `src/activation` | phase 5 — activating function, thresholds per fascicle |
| `src/safety` | charge / charge-density checks |
| `src/utils` | shared helpers (plotting, I/O, unit conversion) |
| `config/` | parameter files |
| `scripts/` | runnable pipelines (`run_*.m`) |
| `tests/` | verification tests (`matlab.unittest`) |
| `data/electrode` | provided active-site layout |
| `data/patients` | virtual patients (seed + parameters) |
| `results/` | generated outputs (not versioned) |
| `report/` | final report / slides |

## Running

From MATLAB, in the project root:

```matlab
setup_paths            % adds src/ and config/ to the path
```

Requires MATLAB R2025a with the **Partial Differential Equation Toolbox**.

## 3D geometry preview

From the repository root in MATLAB:

```matlab
setup_paths
[fig, model] = run_3d_quicklook("deflect");
[figStraight, modelStraight] = run_3d_quicklook("penetrate");
```

Rotate or zoom the 3D axes with the figure toolbar. The end-on view overlays the
four mirrored station pairs from `sites.mat`. There are eight axial stations,
with four internal pads and two opposing cuff pads per station. All patches
represent finite active faces (0.2 × 0.2 mm in the supplied file).

- **Inputs** are in `config/quicklook3DParams.m`. The 60 mm length, circular outer
  section and ten straight fascicles are provisional. The illustrative fascicle
  layout was chosen to exercise both interaction scenarios.
- **Nerve and cuff apposition** are in `src/geometry/buildQuicklookNerve.m`. The
  cuff faces are rigid planes. The sampled nerve surface meets each plane over
  the full pad footprint and returns smoothly to the circle along the axis.
  The flattening depth is prescribed, with no tissue mechanics or conservation
  of cross-sectional area.
- **Strip and pads** are in `src/electrode/planQuicklookStrip.m`. The
  `deflect` case searches a family of smooth paths between fixed entry and exit
  points, keeping the fascicles fixed. It checks all polyline segments and pad
  projections against the perineurium with a provisional clearance envelope.
  A station with no path in the search family raises an error. The `penetrate`
  case retains the straight path and returns the overlapping fascicle IDs in
  `modelStraight.strips(k).envelopeFascicles`. Overlap is geometric; this mode
  does not model rupture or estimate penetration probability.
- **Visualization** is in `src/utils/plotQuicklookNerve.m`. The sampled nerve
  shell, strip centrelines and active patches can be inspected separately from
  the existing 2D FEM. Mesh generation, field solving and mechanical
  displacement are outside this preview.

The [FINE experiments](https://pubmed.ncbi.nlm.nih.gov/12611367/) provide an
example of a cuff reshaping a nerve. The
[SPINE experiments](https://pubmed.ncbi.nlm.nih.gov/9086385/) report that blunt
elements can separate epineurial tissue while preserving fascicles. A
[penetrating cuff experiment](https://doi.org/10.1186/s42234-025-00193-6)
used guided insertion between fascicles and confirmed interfascicular
placement by micro-CT in a rodent sciatic nerve. These studies use different
devices and do not supply a penetration frequency for this strip in the human
median nerve. The two modes represent geometric scenarios while insertion
mechanics and frequency remain uncalibrated.
