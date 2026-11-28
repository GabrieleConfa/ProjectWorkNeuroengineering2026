# Intraneural stimulation of the human median nerve — FEM project

Neuroengineering course project (Politecnico di Milano, A.Y. 2026/2027) — group work,
tutor Luca Sbarsi (luca.sbarsi@polimi.it).

Does the active-site layout of a new intraneural electrode reach **every fascicle** of the human
median nerve? Answered with a 3D finite element model of the nerve and of the 48 contacts (MATLAB,
mesh by Gmsh, own linear-tetrahedra assembly with the perineurium as a thin resistive layer), an
activating-function criterion for fibre activation, a safety check, and virtual patients. A second
study asks whether the distal contacts can **block** the action potentials the proximal contacts send
towards the hand (MRG fibre model in NEURON, Python notebooks): see `docs/block.md`.

Brief: [docs/project_brief.md](docs/project_brief.md) ·
Architecture: [docs/architecture.md](docs/architecture.md) ·
Work plan: [docs/work_plan.md](docs/work_plan.md) ·
Parameters: [docs/parameters.md](docs/parameters.md) ·
Decisions: [docs/decisions.md](docs/decisions.md) ·
3D model explained: [docs/model_3d.md](docs/model_3d.md) ·
Status and next steps: [docs/design.md](docs/design.md) ·
Meetings: [docs/meetings/](docs/meetings/) ·
Paper notes: [docs/literature/](docs/literature/)

**New to the project?** Read this README, then `docs/model_3d.md` (the 3D model and its code, step by
step) and `docs/block.md` (the block study); the choices and their reasons are in `docs/decisions.md`.
Run `setup_paths; runtests('tests')` to check your installation. The entry points are
`src/main/main_pipeline_3d.m`, `scripts/export_viewer_data.m`, `scripts/export_block_profiles.m` and the
notebooks. Nine functions in `src/` (and `scripts/run_pipeline.m`) are skeletons of the first work plan
that were never implemented; their header says so and points to the code that does the job.

## Structure

| Folder | Content |
|---|---|
| `docs/` | brief, architecture, work plan, parameter table, decision log (`decisions.md`), 3D model walkthrough (`model_3d.md`), block study (`block.md`), status (`design.md`), patient file format (`patient_data.md`), one note per paper (papers themselves are not in the repo: download them from the DOIs in the notes) |
| `src/geometry` | phase 1 — nerve cross-section generator (oval nerve, fascicle layout, cuff pad geometry) |
| `src/patient` | virtual and measured patients (JSON files in `data/patients/personal/`, format in `docs/patient_data.md`) |
| `src/insertion` | Monte Carlo of the lead insertion: where the bent intraneural lead ends up in each station |
| `src/materials` | phase 2 — skeleton only, never implemented: the conductivities are in `config/pw6_3d_params.m` and assigned in `build3DModel.m` |
| `src/electrode`, `src/mesh` | phase 3 — active sites, 3D geometry and mesh (Gmsh mesher `meshPersonalDevice.m` + `mesh_personal.py`, default; PDE Toolbox mesher `buildNerveGeometry3D.m`, legacy) |
| `src/fem` | phase 4 — P1 assembly, perineurium thin layer, interpolation, grid location |
| `src/activation` | phase 5 — skeletons only, never implemented: activation and thresholds were done with the MRG fibre model (`src/block/`, `notebooks/`) |
| `src/safety` | charge / charge-density and field-exposure checks |
| `src/block` | conduction-block study: FEM potential along the fibres (MATLAB) and the MRG fibre model, batch runner, volley analysis and anode design (Python + NEURON) |
| `src/utils` | shared helpers |
| `src/main` | the runnable pipelines: `main_pipeline.m` (2D), `main_pipeline_3d.m` (3D simulation, calls `build3DModel.m` and `solve3DContacts.m`) |
| `config/` | parameter files: `pw6_3d_params.m` (3D model), `block_params.m` (block study), `quicklook3DParams.m` (geometry of the Gmsh model), ... — the only place where numbers are set |
| `scripts/` | helper scripts: viewer export (`export_viewer_data.m`), fibre profiles for NEURON (`export_block_profiles.m`), convergence, field safety, insertion study, MATLAB 3D view; `notebook_generators/` writes the block notebooks |
| `notebooks/` | Jupyter notebooks of the block study (Python + NEURON); each one starts with what it asks and how to run it |
| `tests/` | verification tests (`matlab.unittest`): analytical FEM cases, geometry, patients, insertion, safety |
| `data/electrode` | provided active-site layout |
| `data/patients` | virtual patients (seed + parameters) and personal patients (JSON) |
| `results/` | generated outputs; only the viewer data (`results/viewer/`) and the 3D figures (`results/main3d/*.png`) are versioned |
| `report/` | final report / slides |
| `viewer/` | interactive 3D viewer in the browser (`index.html`, three.js vendored in `viewer/vendor/`) |

## Running

Requires MATLAB R2025a with the **Partial Differential Equation Toolbox** (check with `ver('pde')`).
Every command below is typed in MATLAB, in the project root.

```matlab
setup_paths                 % adds src/, config/ and scripts/ to the path (once per session)
runtests('tests')           % verification tests: all must pass before pushing
```

### Just looking at the results: the 3D viewer (no MATLAB needed)

The data of the last simulation (patient UA01, oval nerve, Gmsh model with the bent leads, all 48 contacts
solved, lateral surface grounded and end planes insulated) and 13 exported NEURON runs of the block study are in
the repository, so after cloning you can **double-click `viewer/index.html`**. It opens in Chrome, Edge or
Firefox, from disk, also offline.

| Do this | To |
|---|---|
| drag with the left / right button, scroll | rotate / pan / zoom the nerve |
| hover a gold contact | see where it is (fascicle or epineurium) and its access resistance |
| click a contact, or press **Space** | stimulate it: its field appears as a glowing heatmap (Space again = off) |
| **← →** | previous / next solved contact |
| **1 2 3 4 5** or the buttons | show V, \|E\|, \|J\|, **AF+** (activating function, depolarisation, the default) or **AF−** (hyperpolarisation, patterns only) |
| *Pattern* button, or **Shift+click** contacts | stimulate **several contacts at once** with their own currents (see *Patterns* below) |
| **C**, then **↑ ↓**, **F** | cut the nerve (none → transverse → lengthwise), move the cut, flip the side removed |
| station buttons −4 … +4 | jump to the cross-section through that station |
| **P**, **M**, **L**, **R**, **?** | pulse animation, colour map, full 70 mm, reset view, list of keys |

What the four quantities show (all for a monopolar stimulus of 1 µA; the "Current" box rescales
the colour bar, the model being linear):

- **V**, potential [mV]: highest at the contact, but the whole nerve sits at a few mV because the
  insulating sleeve forces the current along the nerve. The potential of the contact (= its access
  resistance × I, shown when hovering it) depends on where it sits: shank contacts 9–18 kΩ, lowest
  when the contact lies inside a large fascicle (e.g. AS 4, half in F1: 11.9 kΩ), cuff pads 27–39 kΩ
  (one face only, and the current must run along the insulated nerve).
- **|E|**, field [V/m], and **|J|**, current density [µA/mm²], both on a log scale: |J| shows where
  the current flows. Current injected inside a fascicle is kept there by the perineurium and runs
  along the fascicle.
- **AF**, activating function [V/m²]: −∂²V/∂z² over one internode (1 mm), i.e. for a cathodic pulse,
  shown only inside the fascicles and only where it depolarises. **This is what decides whether a
  fibre fires**, so it is the view that answers the project question: a shank contact drives the
  fascicles next to it far more than the others, the cuff pads spread over the nerve surface
  (`export_viewer_data` prints, for every contact, the activating function reached by 10 % of the fibres
  of each fascicle).

Colours are scaled per contact (each contact's own range): to compare contacts, read the numbers
on the colour bar. A link can open a prepared view, e.g.
`viewer/index.html?site=2&heat=1&cut=transverse&z=-10.3&flip=1` (options: `site`, `heat=1`,
`q=V|E|J|A|N`, `cut=transverse|longitudinal`, `z`, `angle`, `offset`, `flip=1`, `full=1`,
`cmap=inferno`, `patient`, `pattern=id:µA,...`). Details: [docs/model_3d.md](docs/model_3d.md) §8.2.

#### Patterns: several contacts at once (2026-11-26)

Press **Pattern** in the *Stimulation* panel. A table appears with one row per station (−4 … +4) and one box per
contact (the grey number is the contact id). Type a current in µA in any number of boxes: **negative = cathodic**,
positive = anodic, empty = off (cathodes turn cyan, anodes red). The fields are summed in the browser as soon as
you leave a box (or press *Apply*); **Shift+click** a contact in the 3D view cycles it through −1 → +1 → off µA.
*Invert* swaps cathodes and anodes, *Clear* switches everything off. The model is linear, so a pattern is the sum
of the single-contact potentials weighted by their currents (exported per contact as a 16-bit potential):

- **V** is coloured by its magnitude **|V|** (cathodes give V < 0: a signed min→max scale would put the
  contact at the bottom and light up the whole nerve); the signed range is printed under the table;
- **V** and **AF+ / AF−** are computed from that sum: AF = ∂²V/∂z² over one internode with the actual polarities,
  inside the fascicles; AF+ shows where fibres are depolarised (pushed to fire), AF− where they are
  hyperpolarised (e.g. under an anode, the block region, and the virtual anodes beside a cathode);
- **|E|** is a finite difference on the 50/100 µm grid (approximate next to the perineurium, where V jumps);
  **|J|** is available for single contacts only;
- the colour bar shows absolute values for the currents you typed (the *Current* box is not used).

A link can open a pattern: `viewer/index.html?pattern=2:-10,44:20&heat=1` (contact:µA, comma-separated).
With the Gmsh model (D-14) the viewer also draws the irregular fascicles displaced by the leads, the bent PDMS lead
of each station (*leads* checkbox) and the PDMS sleeve (*cuff / sleeve*); the cut faces show the exact outlines.
Data exported before 2026-11-26 still open (single contacts, circles), but patterns need a new `export_viewer_data`.

#### Volley playback: the NEURON runs in slow motion (2026-11-27)

When runs of `notebooks/block_volley.ipynb` have been exported (last cell of the notebook, or
`volley.export_viewer_volley` in `src/block/volley.py`), a **Volley · NEURON playback** section appears in the panel.
Nothing is simulated in the browser: it replays the action-potential times recorded by NEURON.

- **Run**: one simulated scenario, e.g. the proximal pulse alone, proximal + anode, or the anode alone;
- **Play** (or <kbd>V</kbd>): every fibre that fires is drawn as a glowing string in its fascicle's colour, and each AP
  as a bump travelling along it (**Style** *wave · string*), or as an oscillation under a moving envelope
  (*wave · ripple*, like a wave packet); *nodes · dots* shows instead the nodes of Ranvier lighting up. The bump
  follows the AP times recorded by NEURON (its length on screen is conduction velocity × AP duration, ~3–4 mm);
  its **Height** (outwards from the nerve axis) is only a drawing, not a physical quantity, and grows with the length
  shown so the bump stays visible at full length. **Colour**: neon cyan (default), magenta, lime or yellow, with a
  white-hot crest and a halo of fixed size on screen (visible at any zoom), or the fascicle colours. **Glow** scales
  the brightness (0–3×), **Core** and **Halo** the sizes of the bright centre and of the soft halo in screen pixels
  (0 hides that layer); with many fibres firing, lower Glow or Halo if the crests merge into a white blob. The wave is drawn on
  top of the glass shells; untick *nerve* and *fascicles* (or open the link with `&hide=nerve,fascicles,cuff,leads`)
  to see it best. **Time** scrubs, **Speed** sets how much nerve time passes per second (default 0.2 ms/s: the
  volley crosses the electrode zone in a few seconds);
- *all fibres* draws every fibre line faintly; *target only* hides the other fascicles; turn the *Stimulate* heatmap
  off and switch on *full length* (<kbd>L</kbd>) to follow the volley to the fibre ends (the APs are counted 5 mm
  from each end);
- the chart (bottom left) shows the **volley of each fascicle along the nerve** (fraction of its fibres with an AP
  at each node; dotted line = 10 %) and, below, the stimulator currents in time with a cursor; the panel lists
  how many fibres of each fascicle have delivered an AP to the brain end and to the muscle end so far (orange =
  10 % or more of the fascicle reached the muscle, criterion C1/C2 of D-16).

Links: `viewer/index.html?volley=<name>&vt=0.3&full=1` opens a run at t = 0.3 ms; add `&vplay=1` to start playing,
`&vstyle=string|ripple|dots` to choose the style, `&vcol=cyan|magenta|lime|yellow|fascicle` the colour, `&vglow=1.5&vcore=9&vhalo=30` the glow, `&hide=nerve,fascicles,cuff,leads,contacts,labels` to hide parts.

#### Tutorial: a first look in five steps

1. Open `viewer/index.html`, choose **AS 2** in *Contact* and press **Space**. The fascicles light up
   with the activating function of that contact.
2. Press **C** (transverse cut), then click the station button **−4**: you now see the cross-section
   through the ring of contacts 1–6, with the heatmap painted on the cut face. The cut face is the
   best place to *read* the map, because the 3D glow adds up colour along the line of sight.
3. Move the cut with **↑ ↓** (or the *z* slider) to see how far along the nerve the stimulus reaches.
4. Press **← →** to go through contacts 1–6 and compare which fascicle each one drives.
5. Read the values by matching a colour on the cut face with the colour bar. Type a current in
   *Current* to read the values for that current instead of 1 µA. Hovering a contact shows where it
   is (fascicle or epineurium) and its access resistance.

#### The control panel

| Control | What it does |
|---|---|
| *Contact* | which solved contact is shown |
| *V / \|E\| / \|J\| / AF* | which quantity is shown (see above) |
| *Current* | multiplies the **numbers on the colour bar** only. The model is linear, so the picture is the same at any current |
| *Window lo* | values below it are not coloured (transparent in 3D, tissue colour on the cut face) |
| *Window hi* | values above it all get the top colour |
| *Glow* | opacity of the 3D heatmap; purely visual |
| *Colours* | colour map (turbo or inferno) |
| *pulse train* | animation of a pulse train; purely visual |
| *z*, *angle*, *offset*, *flip side* | position of the transverse or lengthwise cut, and which half is removed |
| *keep contacts* | the contacts are drawn even where the cut removes the nerve |
| *Show* checkboxes | hide/show nerve, fascicles, cuff, contacts, labels; *full length* shows all 70 mm |

**How the window works.** Each contact's field is stored on a 0–1 scale *of its own*.
For V the scale is linear, from the minimum to the maximum. For |E|, |J| and AF it is
**logarithmic over 3 decades**: 1 is the 99.9th percentile of that contact's field, 0 is 1000×
smaller, so **0.33 on the slider = a factor of 10**. The colour bar always shows the window in
physical units: its left end is the value at *lo* and its right end is the value at *hi*.

| Window lo (AF, \|E\|, \|J\|) | Shows values down to | Use it to |
|---|---|---|
| 0.9 | 1/2 of the maximum | see only where the stimulus is strongest (the fibres that would fire first) |
| 0.8 | 1/4 | see the main target fascicle and its closest neighbours |
| **0.5** (AF default) | 1/30 | general view: target and the clearly secondary fascicles |
| 0.6 (\|E\|, \|J\| default) | 1/16 | general view of field and current paths |
| 0.2–0.3 | 1/250–1/125 | see how far the weak spread reaches (the whole nerve lights up) |

Keep *Window hi* at 1 (lowering it saturates the top of the scale). For V the defaults (0.25 / 1)
are fine. *Glow* between 1 and 3: lower it if the 3D view looks foggy, raise it for faint fields.
Once a threshold is known (phase 5), set *Window lo* so that the left end of the colour bar reads
the threshold at your current: everything coloured is then above threshold.

**Reading the colours correctly.** Colours are relative: red means "close to the strongest value of
*this* contact", not "the fibres fire". Whether fibres fire depends on the absolute value (colour bar
× current) compared with a threshold. It is expected that the fascicle nearest a contact is the
most stimulated, even through the epineurium. The epineurium is a conductor (0.083 S/m). The
resistive barrier is the perineurium (0.00088 S/m, 14 µm), and it surrounds *every* fascicle. It
lowers the current that enters all of them and so raises the threshold, but it does not change
which fascicle is reached first: the potential falls with the distance from the contact.

### Running the 3D simulation

```matlab
setup_paths; main_pipeline_3d       % or open src/main/main_pipeline_3d.m and press Run
```

The default (Gmsh) mesher calls Python: it needs `python -m pip install gmsh numpy scipy`, or the
environment variable `NERVE_GMSH_PYTHON` pointing to a Python that has them.

What it does, in order (explained step by step in [docs/model_3d.md](docs/model_3d.md)):

1. generates the virtual patient (oval nerve 4.35 × 2.69 mm, roundness 0.62; 10 fascicles with
   irregular outlines, mid-upper arm, from the literature) and saves it in `data/patients/<id>.mat`;
2. simulates one random insertion of the PDMS lead per station (it bends around or pierces the
   fascicles, `src/insertion/`) and places the four shank pads along it; the cuff pads are rigid and
   flat, tangent to the oval and pressed in until the whole face touches the nerve; the nerve
   (70 mm) is wrapped in an insulating PDMS sleeve (50 mm);
3. meshes it with linear tetrahedra through Gmsh (minutes, only the first time: the mesh is saved in
   `results/mesh/` and reused; never use `p.mesh3d.cacheMode = "regenerate"` for production runs,
   it invalidates every solved contact);
4. assembles the FE problem (perineurium as a thin resistive layer) and solves it once per contact
   with 1 A (pcg, about 30 s per contact; each contact is saved as soon as it is solved and a re-run
   resumes from what is on disk);
5. saves `results/main3d/<id>_3d.mat` (mesh + potentials) and the figures `results/main3d/*.png`.

Every choice is a parameter in **`config/pw6_3d_params.m`** (never in the code); the ones you will
change most:

| Parameter | Default | Meaning |
|---|---|---|
| `p.patient.id`, `p.patient.seed` | `"UA01"`, `1` | which virtual patient (same seed = same patient) |
| `p.geometry.generator.nerveShape`, `.nerveRoundness`, `.nerveRotation_deg` | `"ellipse"`, `0.62`, `0` | nerve section (`"circle"` = the first, round model) and orientation of the oval in the electrode frame |
| `p.mesh3d.mesher` | `"gmsh"` | `"gmsh"` (bent leads, irregular fascicles, D-14) or `"pdetoolbox"` (straight shanks, circles, the earlier model) |
| `p.solve3d.sections` | `[]` | stations whose contacts are solved; `[]` = all 48 contacts |
| `p.domain.endBC` | `"insulated"` | end planes `"insulated"` (D-17) or `"ground"`; the lateral surface is grounded |
| `p.mesh3d.h*_m` | 5 µm (pad edges) … 1.2 mm | element sizes |

Then, to see the new results in the viewer:

```matlab
setup_paths; export_viewer_data     % a few seconds per contact
```

It rewrites `results/viewer/<id>/` (one `AS<k>.js` per solved contact) and prints, for each contact,
the activating function reached by 10% of the fibres of each fascicle. Reload the viewer page.
Commit the new `results/viewer/` files if the group should see them (about 7 MB per contact).

Other entry points: `src/main/main_pipeline.m` (2D model), `scripts/view_model_3d.m` (static 3D view
in MATLAB), `scripts/run_convergence_3d.m` (mesh convergence), `scripts/run_field_safety.m` (field
exposure), `scripts/run_insertion_3d.m` (lead-insertion statistics).

### The conduction-block study (Python + NEURON)

Question: the proximal contacts excite a fascicle; can the distal contacts stop the action potentials
that run towards the hand? Physics, pipeline, every run and its result: [docs/block.md](docs/block.md);
choices: D-12, D-15 … D-20 in [docs/decisions.md](docs/decisions.md).

Requirements: NEURON (Windows installer, `C:\nrn`; `pip install neuron` has no Windows wheels), the MRG
membrane mechanism compiled once in `results/mrg/model/` (`nrnivmodl`, see `docs/block.md` §4), and
`h5py`, `numpy`, `scipy`, `pandas`, `matplotlib`. Then:

1. in MATLAB, after every new 3D solution: `setup_paths; export_block_profiles` (FEM potential along
   ~450 fibre lines, `results/block/<id>_fibre_profiles.h5`, and a JSON copy of `config/block_params.m`);
2. open a notebook in `notebooks/` and *Run All*. `block_volley.ipynb` is the main one (travelling
   volley, anodal block, whole-nerve check, export to the viewer); `block_anode_design.ipynb` designs
   the minimum anode configuration by optimisation; `block_model_check.ipynb` checks the fibre model
   against textbook DC block; `block_anode_cathode.ipynb` and `block_anode_configurations.ipynb` are the
   earlier rounds. `scripts/run_block_study.py` is the first, command-line version of the study.
