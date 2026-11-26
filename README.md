# Personal median nerve model · MATLAB

Open **start_nerve.m** in MATLAB R2025a and press **Run**. The default run regenerates the geometry, figures, connected-section mesh, all 48 electrical solutions and the offline browser views. The Command Window prints the local links.

| To reopen a saved result | MATLAB command | File to double-click after extracting the project |
|---|---|---|
| All three saved pages | `start_nerve("open")` | Open the three files below |
| Stimulation, initially on AS5 | `start_nerve("viewer")` | `viewer/personal.html` |
| Rotate and section the tetrahedral mesh | `start_nerve("mesh")` | `viewer/mesh.html` |
| Review the figure gallery | `start_nerve("figures")` | `viewer/figure.html` |

### Inputs and calculation

- Edit `config/quicklook3DParams.m` for the oval, fascicle coordinates, lead mechanics and cuff closure. Electrical settings are in `config/pw6_3d_params.m` and the supplied active sites are in `data/electrode/sites.mat`.
- Each supplied cuff azimuth identifies a point on the native oval. The flat face is tangent to the oval at that point before closure and advances along its local normal. The contact figure shows two pads on different parts of the oval. The 3 µm extra closure per pad is an assumption. The surrounding insulating sleeve uses a separate assumed circular outer envelope.
- Ten representative fascicles occupy 5.01 of 9.20 mm². The largest is 1.49 mm². The ellipse represents group measurements. Individual centres, other areas and functional assignments are unknown. The provided synthetic layout was selected for geometric feasibility of eight insertions. It cannot estimate penetration frequency. [Pušnik et al.](https://doi.org/10.1038/s41598-024-84396-y)
- MATLAB constructs successive nerve contours spanning −30 to +30 mm. The Gmsh helper in `src/mesh` joins corresponding profiles, intersects the device and tissue volumes, and generates conforming tetrahedra. MATLAB assembles and solves the anisotropic conduction problem. To rebuild, install Gmsh, NumPy and SciPy into a Python environment (for example `python3 -m venv .venv` followed on macOS by `.venv/bin/python -m pip install gmsh numpy scipy`). The environment variable `NERVE_GMSH_PYTHON` can specify the Python executable on another computer.
- In the stimulation viewer select a contact and press **Space**. Keys **1–4** select potential, field magnitude, current density and activating function. **P** animates a pulse. The activating function estimates axial electrical drive; the animation does not calculate action potentials. The 60 mm FEM domain includes the ends; the exported heatmap covers the central 25 mm. [Original viewer and solver design](https://github.com/GabrieleConfa/ProjectWorkNeuroengineering2026/blob/9fdf9ea/docs/design.md)
- Helpers live in `scripts/` and `src/`. Generated figures, numerical solutions and browser data live in `results/3d_preview`, `results/main3d` and `results/viewer`. The HTML files use these relative paths, so keep the folders together.
- `start_nerve("convergence")` repeats the mesh study (four meshes, pads AS1, AS5, AS17) and writes tables and figures to `results/convergence`; the gallery shows them. The study reuses the production mesher and solver and does not overwrite the 48 saved fields.
- Older outputs (2D FEM, local slab, reference-patient figures, single-fibre MRG runs) were moved out of the project to `../archivio_risultati_vecchi/results/`. To rerun the single-fibre experiment below, move `mrg` back into `results/`.
- `start_nerve("insertion")` runs the scenario study with random irregular fascicles and Monte Carlo lead insertions (`scripts/run_insertion_study.m`, parameters in `config/insertionParams.m`, algorithms in `docs/algorithms_and_plan.md`). It is not yet part of the 3D mesh and viewer. The mesh uses 5 µm elements on all pad edges (`p.mesh3d.hEdge_m`), tested by `scripts/run_edge_refinement_test.m`.
- Geometry mode: `p.geometry.mode` in `config/quicklook3DParams.m`. `"insertion"` (default) uses random irregular fascicles and one random PDMS lead insertion per station, from the entry wall to the opposite wall; the mesh, the solver, the figures and the viewer read the same scenario (`src/insertion/insertionScenario.m`), and the mesh contains PDMS lead bodies and an oval PDMS sleeve. `"legacy"` restores the circular fascicles and the previously approved lead routes.
- Pending measurements are listed in [the meeting notes](docs/meetings/2026-11-23_domande_prof.md).

<details>
<summary>Group reference workflows and development notes</summary>

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
3D model explained: [docs/model_3d.md](docs/model_3d.md) ·
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
| `viewer/` | interactive 3D viewer in the browser (`index.html`, three.js vendored in `viewer/vendor/`) |

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
- **Nerve and cuff apposition** are in `src/geometry/cuffContactPlane.m` and
  `src/geometry/cuffDeformedRadius.m`. Each rigid plane is fitted to its finite
  transverse pad width from `sites.mat`. The smallest central indentation that
  puts the original circular or oval contour in contact with the whole pad is
  computed from the contour and the pad edges, not entered as a fixed depth.
  Only the nerve surface facing that pad is flattened. Smooth transverse and
  axial returns outside the pad have editable, provisional widths of 0.12 and
  0.40 mm. This geometric construction does not solve tissue mechanics,
  conserve cross-sectional area or predict movement of internal fascicles.
- **Strip and pads** are in `src/electrode/planQuicklookStrip.m`. The
  `deflect` case searches a family of smooth paths between fixed entry and exit
  points, keeping the fascicles fixed. It checks all polyline segments and pad
  projections against the perineurium with a provisional clearance envelope.
  A station with no path in the search family raises an error. The `penetrate`
  case retains the straight path and returns the overlapping fascicle IDs in
  `modelStraight.strips(k).envelopeFascicles`. Overlap is geometric; this mode
  does not model rupture or estimate penetration probability.
- Each station has one continuous lead entering from the negative side of the
  local insertion axis and ending beyond the fourth active face. The horizontal
  station enters from the left. Other stations rotate with the supplied pad
  coordinates. The nominal tip extension is provisionally 50 µm
  (`p.electrode.tipOverhang_m`). The FEM represents the conductive contacts
  and the three tissues. The carrier path is a geometric overlay.
- **Visualization** is in `src/utils/plotQuicklookNerve.m`,
  `src/utils/plotStationAtlas.m` and `src/utils/plotCuffApposition.m`. The 3D
  nerve, all eight current transverse stations, the strip route, and a
  micrometre-scale view of local cuff contact can be inspected separately from
  the existing 2D FEM. Field solving and mechanical displacement are outside
  this preview.

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

The experimental basis and limits of the insertion scenarios are recorded in
[`docs/insertion_evidence.md`](docs/insertion_evidence.md). A first local
three-dimensional tetrahedral mesh of one flattened cuff station can be made
with `buildQuicklookMesh(p, cuffSites)`. It is a short uniform slab through a
contact. The separate full-length personal case is described below.

The insertion screen now reports a straight-path encounter angle and compares
an ideal weak-axis buckling load with a probe-specific perineurium insertion
force. Its dimensions, unsupported length and material proxy are editable in
`quicklook3DParams.m`; changing the strip width also changes the geometric
collision envelope. The screen cannot predict a penetration percentage or
confirm rupture. The source and assumptions are listed in
`docs/insertion_evidence.md`.

The current 3D anatomy is entered in `config/quicklook3DParams.m` as the
nerve's two transverse semi-axes, rotation, and one `[x y semi-axis-x semi-axis-y rotation-deg]`
row per fascicle, all in SI units. The 3D view accepts a different fascicle
count and an oval nerve outline. `validateQuicklookAnatomy` rejects invalid
sizes, clear boundary violations and overlapping fascicles. The local slab
mesh uses the same locally indented contour and accepts editable oval nerve
and fascicle axes. It represents only the 0.20 mm axial contact slab, so it
does not contain the axial contact shoulders or the full 6 cm nerve. The
historical 2D FEM still uses its earlier circular section.
The original 2D figures are preserved under `results/2d_sites` and describe
the tutor P00 section with seven fascicles. The new 3D figures under
`results/3d_preview` describe the current ten-fascicle quicklook anatomy.
Use `run_complete_preview` to generate the personal geometry and section atlas. The old 2D
solver reports distorted elements and places its cuff sites in muscle in its
own section; its field plots must not be interpreted as results for the new
ten-fascicle model.

The quicklook collision check uses each fascicle's larger semi-axis as a conservative circular envelope. The displayed fascicle is oval when its axes differ. Exact oval contact mechanics remains part of the later calibrated model.

## Local mesh and electrical verification

```matlab
setup_paths
[meshPreview, field] = run_local_mesh_stimulation();
```

This produces an inspectable MATLAB mesh figure, a cross-sectional potential
map, and the numerical field under `results/local_slab/`. The source is a
finite cuff pad held equipotential; the opposing pad is the return. The field
is scaled to 1 µA and uses the project's anisotropic endoneurium, perineurium,
and epineurium conductivities. It is a **short-slab solver check**, not a
prediction of stimulation of the full nerve or of action-potential propagation.
The pad boundary is represented by mesh nodes, and the present mesh is not
fully convergence-validated. A second mesh with 188,827 tetrahedra changed
the local two-terminal conductance by about 0.06%; this checks one quantity
for one station only. All eight stations have been run separately with the
same local setup. The full 3D pipeline uses the MATLAB R2025a geometry API.

### Interactive 3D viewer

After `main_pipeline_3d` has solved some contacts:

```matlab
setup_paths; export_viewer_data     % writes results/viewer/<patient>/ (geometry + one file per contact)
```

then open `viewer/index.html` in Chrome, Edge or Firefox (double-click, no server needed). Rotate the
nerve with the mouse, press Space (or click a contact) to show its field, cut the nerve
transversally or lengthwise. Keys and options: `docs/model_3d.md` §8.2.

### Personal 60 mm case

From the project root in MATLAB R2025a, run `setup_paths; run_personal_3d`.
The original reference patient is UA01. The personal case is Cuff60 and uses the
48 contact IDs from the supplied device. The eight section drawings and the
finite-element solver use the same routed internal pad centres and orientations.
The personal mesh uses Gmsh 4.15.2 with NumPy and SciPy. Install those packages
in a Python environment and set `NERVE_GMSH_PYTHON` to its Python executable.
The default is `~/.cache/nerve-gmsh/bin/python`. MATLAB fits the surface
nodes to the local cuff displacement, checks signed tetrahedral volumes and
measures the triangles shared by each active pad and the epineurium before
solving. Full apposition is an imposed geometric condition. The cuff pressure
and tissue response still require measurements.
Open `viewer/personal.html` for the stimulation view. Field sampling for
display uses 100 µm in x and y and 200 µm in z. The numerical solution retains
the tetrahedral mesh resolution. The monopolar return is the grounded exterior.
For dimensions and functional assignments, set `p.geometry.functionLabels`
alongside the fascicle rows and run `plotFascicleMap(p)`.

### Single-fibre propagation experiment

The model is the [published MRG implementation on ModelDB](https://modeldb.science/3810),
revision b15daeda6ea368df0198d257f339cf5e840cf31c. The example uses one
10 µm motor fibre through one fascicle centre. The fixed endpoints and the choice
of a single central fibre can affect initiation and recruitment. Recordings at
four nodes and a zero-current control distinguish a voltage response from a
local extracellular field. The supplied pulse is an exploratory simulation input.

In a terminal, from the project root, prepare the Python environment once:

```sh
mkdir -p results/mrg                            # Keep downloaded and compiled files with local results.
cd results/mrg                                 # NEURON discovers compiled mechanisms in this folder.
git clone https://github.com/ModelDBRepository/3810.git model  # Obtain the published model.
git -C model checkout b15daeda6ea368df0198d257f339cf5e840cf31c # Pin the morphology and ion mechanisms.
python3 -m venv env                            # Isolate the propagation tools.
env/bin/pip install neuron numpy matplotlib    # Install the numerical and plotting dependencies.
env/bin/nrnivmodl model                        # Compile the published AXNODE mechanism.
```

From the MATLAB project root, export a solved axial profile:

```matlab
setup_paths % Locate the shared interpolation routine.
export_axon_profile('results/main3d/Cuff60_3d.mat',1,1,'results/mrg/Cuff60_AS1_F1.csv') % Select AS1 and the F1 centreline.
```

From the terminal in `results/mrg`, run
`env/bin/python ../../scripts/run_mrg_fibre.py Cuff60_AS1_F1.csv --current-uA 100 --phase-ms 0.05`.
This writes the voltage traces and the plot with the pulse settings in the filenames.

</details>
