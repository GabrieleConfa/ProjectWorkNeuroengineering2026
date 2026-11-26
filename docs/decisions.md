# Decision log

One entry per modelling choice: what was decided, the alternatives, the source, and the
consequence for the results (it feeds the "limitations" part of the discussion).

Format:

```
## D-NN — <short title>  (YYYY-MM-DD, phase N)
Decision: ...
Alternatives: ...
Why / source: ...
Consequence / limitation: ...
```

## Open decisions (to settle before the corresponding phase)

The decisions needed for sprint 1 are laid out with options and a proposal (codes D-a … D-h) in
`docs/meetings/2026-11-10_kickoff.md` §3; once agreed, they get a D-NN entry below. Facts already
checked (2026-11-08, `docs/architecture.md` §5): `femodel` accepts only isotropic σ, `createpde`
accepts a diagonal tensor in 3D; `extrude`, `addCell` and `addVoid` work for the nerve + AS
geometry. The 3D solve is required (stacked 2D solutions are not equivalent, architecture §2).

- **Implantation level** (phase 1): upper arm (≈ 12 fascicles, some large) vs proximal forearm
  (≈ 22 small fascicles). Delgado-Martínez 2016 recommends either for intraneural electrodes.
- **How the cross-section is generated** (phase 1): parametric random generator from the
  Delgado-Martínez statistics vs. segmentation of a histological image.
- **Conductivity set** (phase 2): Zelechowski 2020 vs Horn 2023 vs a mix; which to use as
  baseline and which in the sensitivity analysis.
- **Perineurium representation** (phases 2–3): explicit thin layer in the mesh vs thickened layer
  with rescaled conductivity (same ρ·t) vs thin-layer boundary condition (if the solver allows it).
- **Outer domain and boundary conditions** (phase 3): size, ground vs insulated outer surface,
  location of the return electrode.
- **Electrode geometry in the mesh** (phase 3): explicit shaft volume vs AS-only boundary patches;
  whether PDE Toolbox can build the geometry or an external mesher is needed.
- **Activation criterion and threshold** (phase 5): activating-function form and the threshold
  value per fibre diameter.
- **Number of virtual patients and what varies between them** (phase 6).

## Decisions taken

## D-01 — Location and loader of the provided electrode layout  (2026-11-08, phase 3)
Decision: the tutor's `sites.mat` lives in `data/electrode/sites.mat` (moved from a zip-extracted
folder `src/MatlabCode/sites.mat/`, content unchanged, MD5 in `data/electrode/README.md`). The
project code reads it only through `src/electrode/loadActiveSites.m`, which converts to SI units
and computes the cuff radius as R + thickness/2. The tutor's `load_active_sites.m` is kept
untouched as reference (it points to a non-existent `active_sites.mat`).
Alternatives: fix the path inside the tutor's script and use it as is.
Why / source: project layout (input data in `data/`, SI units inside the code); one loader means
one place to change if the tutor sends a new layout.
Consequence / limitation: none on the results; `tests/testLoadActiveSites.m` fails if a new data
file breaks the stated layout (48 AS, mirror symmetry, radii, orientations).

## D-02 — 3D domain: nerve 70 mm, cuff 50 mm centred on the contacts  (2026-11-13, phase 3)
Decision: nerve and surrounding muscle 70 mm long (z = −35 … +35 mm), muscle up to 3 R as in 2D;
cuff 50 mm long, centred on the device (the contacts span |z| ≤ 10.4 mm, so ≈ 14.6 mm from the last
station to each cuff end and 10 mm of bare nerve beyond). Lateral surface of the muscle grounded.
All values in `config/pw6_3d_params.m` (`p.domain`, `p.cuff`).
Alternatives: 40 mm domain (previous proposal); cuff offset from the contacts.
Why / source: group choice (lengths); centring is the natural assumption for a symmetric device.
Consequence / limitation: the end planes are only 10 mm beyond the cuff; their condition (D-04)
must be shown not to change the potential at the stations.

## D-03 — 3D geometry by layered extrusion of one cross-section (not slices + interpolation)  (2026-11-13, phases 1–3)
Decision: one 2D section of the tissues (muscle, cuff wall, nerve, fascicles) is extruded along z
in 5 layers (`extrude` with a vector of heights; planes at the cuff ends and at the fine-mesh zone).
The cuff-wall ring is deleted inside the cuff (void) and merged into the muscle outside. The 48
contact boxes are then inserted in 3D with a Boolean union that keeps the existing boundaries
(`union(..., KeepBoundaries=[true true])`, R2025a); the pieces of a box that crosses a tissue
boundary are merged into one contact cell. Code: `src/mesh/buildNerveGeometry3D.m`.
(First version, abandoned the same day: the contact footprints drawn in the 2D section and a layer
per station, i.e. 21 layers and ≈ 1000 pairwise `mergeCells` — still running after 40 minutes.)
The nerve was circular until 2026-11-22 (the cuff rigid and tight, so the nerve took the round
shape of the cuff); since then it is an ellipse, see D-10. The nerve still conforms to the electrode,
not vice versa (group, 2026-11-13; D-11).
Alternatives: slices at several z joined by interpolation (lofting).
Why / source: (i) no data on how a given fascicle changes along 7 cm — interpolating would invent
it, and fascicles that split or merge cannot be interpolated; (ii) the PDE Toolbox has no lofting
(an external CAD tool / mesher would be needed); (iii) with extrusion a fibre is a straight line
parallel to z, which is what the activating function needs (phase 5); (iv) continuity with the 2D
section.
Consequence / limitation: fascicles do not change along z (real fascicles split and merge every few
mm; Upadhye 2022 on the vagus: every ≈ 560 µm). The variability is put between patients instead.

## D-04 — End planes of the domain grounded (to compare with insulated)  (2026-11-13, phase 3)
Decision: z = ±35 mm grounded (`p.domain.endBC = "ground"`), like the lateral surface.
Alternatives: insulated ends (`"insulated"`).
Why / source: grounded ends mimic a nerve that continues towards a distant ground; insulated ends
mimic a truncated nerve. Neither is exact for a finite model.
Consequence / limitation: TO DO — run both and report the change of V and of the activating
function at the stations; if it is not negligible, lengthen the domain.

## D-05 — Cuff: rigid, tight, insulating wall modelled as a void  (2026-11-13, phase 3)
Decision: cuff inner surface = nerve surface (no saline gap); inside the cuff the wall (from the
nerve outline to the outline grown by the wall thickness) is deleted from the geometry, so the nerve
surface under the cuff is an insulated boundary. Outside the cuff the same ring is muscle. Wall
thickness 0.5 mm `[assumed]`. Since D-10 the cuff is shaped on the native oval: inner surface = the
nerve ellipse a × b, outer surface = the ellipse (a + wall) × (b + wall) (not an exact parallel
curve; the wall is a void, so its outer shape only sets where the muscle starts).
Alternatives: a meshed low-conductivity wall (silicone); a thin saline film between cuff and nerve.
Why / source: tight fit = group choice; a void is the exact limit of an insulator and costs no
elements.
Consequence / limitation: the wall thickness is not known (question added to
`docs/project_brief.md`); it only sets where the muscle starts outside the cuff, so a small effect
is expected — to verify.

## D-06 — Perineurium as a thin layer (contact impedance), own FE assembly  (2026-11-13, phases 2–4)
Decision: the perineurium is not meshed. Each fascicle boundary (circle of radius rEndo) is a surface
across which V jumps, with J = (V_in − V_out)/(ρ·t), ρ·t = t/σ_peri = 14 µm / 0.00088 S/m
= 159 Ω·cm². Implemented by duplicating the interface nodes and coupling the two sides with the
interface mass matrix (`src/fem/thinLayerInterface.m`). The stiffness matrix is assembled by the
project's own P1 code (`src/fem/assembleConductionP1.m`); the PDE Toolbox only builds the geometry
and the mesh. Verified against exact solutions in `tests/testThinLayerFEM.m`.
Alternatives: meshed perineurium (40 µm, σ scaled to keep ρ·t): tested 2026-11-13, **27.7 million
tetrahedra / 4.6 million nodes** for the 70 mm model — not solvable. Thicker layers would eat most
of the small fascicles.
Why / source: Pelot 2019 (a contact-impedance boundary gives thresholds similar to a meshed
perineurium); ASCENT (Davis 2023) uses it. The PDE Toolbox has no interface condition, hence the
own assembly.
Consequence / limitation: linear elements only (quadratic ones would need a new assembler); the
layer has no thickness, so the endoneurium keeps its full anatomical radius.

## D-07 — Contacts: active = equipotential conductor with a volume source; idle insulating; no shank  (2026-11-13, phase 3)
Decision: every contact is a cell of the mesh (50 × 200 × 200 µm, from its footprint and its station
layer). The active contact gets σ = 1e3 S/m and a uniform source I/volume: it carries exactly I and
is equipotential (the pipeline prints the spread of V on it). Idle contacts: σ = 1e-7 S/m
(insulating, as in 2D) or, by parameter, floating metal. Cuff contacts are flat (rigid electrode)
and pressed into the nerve so that the whole face touches it (until 2026-11-22: 10 µm at the centre
on the circle, where a flat 200 µm face departs from it by 2.9 µm at its edges; now placed on the
oval as in D-11). No shank body (group: contacts only for now).
Alternatives: contact faces as boundaries with an imposed current density (not equipotential); the
real gold conductivity (ill-conditioned matrix).
Why / source: same physics as the 2D model (equipotential contact, idle contacts insulating);
1e3 S/m is ≈ 10⁴ times the tissue, enough for equipotentiality without spoiling the conditioning.
Consequence / limitation: the shank (insulating carrier) is missing, so current leaves the shank
contacts in all directions; to add when its size is known.

## D-08 — 3D mesh: linear tetrahedra, fine zone |z| < 12 mm  (2026-11-13, phase 3)
Decision: element size 30 µm on the contacts, 100 µm on the fascicle boundaries and 150 µm on the
nerve inside the fine zone, 500 µm on the nerve outside it, up to 1.2 mm in the muscle (`p.mesh3d`).
Alternatives: uniform mesh (too large).
Why / source: first guess sized on the mesh budget.
Consequence / limitation: TO DO — convergence study (V at the contacts, activating function along
one fibre) before using any threshold.

## D-09 — Upper-arm cross-section generator (10 fascicles, literature-based, not symmetric)  (2026-11-13, phase 1)
Decision: `makePatient` draws a seeded random section (`src/geometry/generateFascicleLayout.m`):
nerve area 9.20 mm², 10 fascicles, total fascicle area 5.01 mm², one large fascicle of 1.49 mm²
(Pusnik 2025, as in PW6); the other nine with log-normal diameters, CV 0.5 (Brill & Tyler 2017:
0.6 ± 0.3 mm in the upper arm), scaled to the total area; positions uniform random, then relaxed
apart (Delgado-Martínez 2016: fascicles spread homogeneously regardless of size); 2 ± 1 inner
fascicles (Brill & Tyler 2017: ≈ 0.58 (N − 6.08)). Gaps of 25 µm between fascicles and 40 µm from
the nerve surface `[assumed, mesh]`. Shank contacts may cross a fascicle boundary ("pierce": the
rigid shank goes through the tissue and cuts the perineurium); only nearly tangent configurations
are redrawn (mesh quality). Patient UA01 = seed 1, saved in `data/patients/UA01.mat`.
Alternatives: "pushAside" (no contact is ever cut by a perineurium): with z-invariant fascicles all
16 shank footprints constrain the same section, and the only layouts left put the big fascicle on
the axis and the others centred on the contacts — electrode-shaped, nearly symmetric, and biased
towards "every contact is in a fascicle". Kept as an option for comparison. The tutor's P00 is
still available (seed = NaN).
Why / source: group request (upper arm, 10 fascicles, no symmetry, biological variability).
Consequence / limitation: circular fascicles (Brill & Tyler: 97% have circularity > 0.75); the
spread of the sizes comes from one study of 4 nerves; whether the shank pierces or pushes the
fascicles aside is still a question for the tutor.
**To revise (tutor, 2026-11-23):** the probe pierces each fascicle with a probability that grows with
its radius and bends around the others; the "pierce always" rule is to be replaced (`design.md` §5.1).

## D-10 — Elliptical nerve section (roundness 0.62), cuff shaped on it  (2026-11-22, phases 1, 3)
Decision: the nerve section is an ellipse with the literature area (9.20 mm²) and minor/major axis
ratio 0.62: semi-axes R/√0.62 × R·√0.62 = 2.173 × 1.347 mm (R = 1.711 mm, the old equivalent
circle). Major axis along x of the electrode frame (`p.geometry.generator.nerveRotation_deg = 0`
`[assumed]`). The fascicle generator (D-09) packs the same 10 fascicles inside the oval, with the
same gaps, measured as distance from the elliptical boundary; patient UA01 (seed 1) therefore has a
new layout. The rigid cuff is tight on the native oval (D-05). `nerveShape = "circle"` gives back
the round nerve. Code: `src/geometry/generateFascicleLayout.m`, `ellipseQ.m`, `nerveOutline.m`,
`src/mesh/buildNerveGeometry3D.m` (decsg ellipse); the viewer reads `nerveAxes` / `nerveRotation`.
Alternatives: (i) circle (the first 3D model: a rigid round cuff would force the nerve round);
(ii) Marco's "personal" model (branch `cuff-oval-viewer-convergence`): same oval, but one fixed,
hand-written fascicle template and a separate Gmsh pipeline — not taken, because phase 6 needs
seeded virtual patients and the group would need Python + gmsh.
Why / source: group (user, 2026-11-22), after Marco's branch. Roundness 0.62 = Pusnik 2025, Table 1,
histological mean at mid-arm, **as cited on Marco's branch: TO VERIFY against the paper (not in
`Papers/` yet)**. An ellipse with b/a = k also has ImageJ roundness 4A/(π·major²) = k, so the value
does not depend on which of the two definitions the paper uses.
Consequence / limitation: the orientation of the device on the nerve is unknown (question for the
tutor; it decides which cuff pads face the flat sides and which the curved ends). The fascicles are
still circles and the section is still the same along z (D-03). The muscle now extends to
3 × the major semi-axis (6.5 mm instead of 3R = 5.1 mm).

## D-11 — Rigid flat cuff pads tangent to the oval; the nerve conforms, the electrode does not bend  (2026-11-22, phase 3)
Decision: each cuff pad keeps its angle θ from `sites.mat` as a landmark on the native outline; its
face is tangent to the outline there (normal = ellipse normal, up to 24° from radial at θ = 45° on
this oval) and is pushed along that normal until both side edges of the 200 µm face are
`p.electrode.cuffEdgeIndent_m` = 7 µm `[assumed]` below the native surface. The centre of the face
is then deeper by the local sagitta (1.4 µm on the flat sides, 6.0 µm at the ends of the major
axis): 8.4–13 µm. In the geometry the pad box replaces the tissue in front of its face (Boolean
union), so the nerve surface under the pad is the flat face. Shank contacts stay at the tutor's
positions ("pierce", D-09): no deflection of the shank. Code: `src/electrode/placeCuffContacts.m`,
test `tests/testCuffContacts.m`.
Alternatives: (i) Marco's deformed surface: the tissue is displaced onto the pad with smooth
shoulders (0.4 mm along z, 0.12 mm across), meshed with Gmsh — same contact face, the difference is
a sliver a few µm deep around the pad, at the cost of a second meshing pipeline; (ii) Marco's
shank "deflect" scenario (the shank bends around the fascicles): rejected by the group ("the
electrodes don't bend"); (iii) pads pressed radially (the old circle rule) — on an oval the face
would not be tangent to the nerve away from the axes.
Why / source: group (user, 2026-11-22): rigid electrode, the nerve conforms. Edge depth 7 µm keeps
the old 10 µm at the centre on a circle and avoids slivers of tissue thinner than ~7 µm (mesh).
Consequence / limitation: no tissue mechanics: the indentation is geometric, not a pressure (Marco
lists the closing force / pad gap as a tutor question). The few µm of tissue under the face are
removed rather than compressed.
**To revise (tutor, 2026-11-23):** the shank probe does bend inside the nerve (100 µm probe, exit on
the line of the entry); alternative (ii) above is back on the table (`design.md` §5.1). The rigid,
tangent cuff pads are not affected.

## D-12 — Conduction block: two mirror sites, MRG fibres in NEURON on the existing 3D solution  (2026-11-25, block study)
Decision: the block of the APs travelling towards the amputated limb (`design.md` §5.3) is studied on
the 3D model as it is (UA01, the solution shown in the viewer), with no new geometry or solve:
- **Sites** = the two halves of the provided layout: sections −4…−1 (z < 0) and +1…+4 (z > 0), which
  are mirror images (same angles at −k and +k). z < 0 is called proximal (shoulder) `[assumed]`.
- **Distal contacts by mirror**: V_k'(x, y, z) = V_k(x, y, −z) for the partner k' of a solved contact k
  (same type, radius, angle, opposite z). Exact for this model (same section along z, D-03, and a
  mirror-symmetric layout); only the mesh noise differs.
- **Fibres**: straight lines parallel to z on the viewer grid, one every 100 µm `[assumed]`, at least
  10 µm inside a fascicle circle; lines through a contact are dropped (453 fibres kept, 23 dropped).
  V sampled every 50 µm over |z| ≤ 33 mm, P1 interpolation on the endoneurium side of the thin layer
  (as the viewer). Code: `src/block/sampleFibreProfiles.m`, `scripts/export_block_profiles.m`.
- **Fibre model**: MRG (McIntyre, Richardson & Grill 2002, ModelDB 3810) rebuilt in Python/NEURON
  (`src/block/mrg_fibre.py`) so that length and diameter are parameters: 60 mm long, centred on
  z = 0, one node at the proximal stimulating station `[assumed]`, D = 10 µm for the first runs;
  the FEM potential is applied as `e_extracellular` of every compartment.
- **Scenarios** (`scripts/run_block_study.py`): proximal monopolar cathodic 50 µs pulse (Zelechowski
  2020) at 1.2 × the 10%-of-fibres threshold of the most excitable fascicle `[assumed]`; (a) collision
  = the same pattern repeated on the distal mirror contacts after a delay (0–0.8 ms); (b) anodal block
  = the distal contacts driven with reversed polarity, 1 ms pulse `[assumed]`, 1–16 × the proximal
  current; (c) cathodal block (supervisor, 2026-11-25: "on the other side the same electrodes work
  as cathodes to cancel the signal") = the distal mirror contacts as cathodes at 2–64 × the
  proximal current, onset delay 0–0.6 ms, width 50/200/500 µs `[assumed]`, swept on the target
  fascicle, the best settings re-checked on all fibres. APs are counted 5 mm from each fibre end;
  success = no AP at the distal end and the same APs at the proximal end as the proximal pulse
  alone. Each pattern names its target fascicle; exact thresholds (bisection) are computed only for
  the target fascicle's fibres, the others are simulated once at the stimulation current (user
  2026-11-25, speed: ~15 min → a few minutes; `--all-thresholds` restores the full table).
Alternatives: ring electrodes added to the cuff (rejected: the sites are the existing contacts,
user 2026-11-25); a MATLAB node-only fibre model (CRRSS, Sweeney 1987: simpler and self-contained,
user chose MRG); the activating function alone (it shows where the membrane is depolarised or
hyperpolarised but cannot tell whether an AP is blocked or collides); solving the distal contacts
(≈ 3.5 min each; needed only if the section changes along z); kilohertz (KHFAC) block (later: needs
the frequency-dependent conductivity, `design.md` §5.5).
Why / source: tutor meeting 2026-11-23 and user 2026-11-25 (pattern repeated distally with a delay;
"anodal now, KHFAC later"; MRG in NEURON); MRG is the fibre model of Zelechowski 2020 (same
pipeline, 50 µs pulse, 10% criterion); course slides 03 p. 23: currents > 8 × threshold block
conduction (cathodal block, seen in the model).
Consequence / limitation: one fibre diameter per run; straight, parallel fibres (no fascicle
branching along the arm, so the "same pattern" excites the same fascicles at both sites, which a
real nerve does not guarantee); the mirror trick fails as soon as the section changes along z (then
solve the distal contacts); multi-contact patterns rely on superposition with insulating idle
contacts (D-07, approximate); node alignment is the best case at the proximal station only, so the
distal threshold of the same fibre differs; quasi-static tissue (the 50 µs pulse has content above
the perineurium corner frequency, Horn 2023); no electrode–tissue interface; the safety of the large
anodal currents (Shannon, CIC) is not checked yet.
**Update 2026-11-26 (supervisor):** in the anode–cathode configuration, the active contacts of the
anode (distal) site are the same as those activated on the cathode (proximal) site, i.e. the mirror
partners of the proximal pattern with reversed polarity (AS 2 → AS 44). Spreading the anode over
other contacts (other shank contacts, cuff pads) is therefore excluded; the free variables are the
anodal current and the pulse shape/timing (trapezoidal ramps). Code:
`notebooks/block_anode_cathode.ipynb`, which derives the anode pattern from `PROX_PATTERN`.
**Result and next step (2026-11-26, group):** with the same contact as anode (AS 44) the fibre is
fired by the virtual cathodes from 17.5 µA (1.7 × I_stim), with or without 0.1–1 ms ramps, and no
block is found up to 2 mA (`docs/block.md` §6). As our own extension beyond the supervisor's rule,
to be shown to the tutor, other distal anode configurations will be compared (all 24 distal contacts
as a "wall", the cuff pads only, whole stations, several stations along z). For this the whole
proximal half is solved (`p.solve3d.sections = -4:-1`, contacts 1–24) and the distal half is taken by
mirror; the solve loop now writes each contact to `results/main3d/<id>_contacts/` and resumes from
there.

## D-13 — Merge of Marco's push 274b485: insertion model, convergence method, pipeline split  (2026-11-26, phases 1, 3, 4)
Decision: from Marco's branch `cuff-oval-viewer-convergence` (commit 274b485, 2026-11-26, his
separate Gmsh "personal" model), the parts the user asked for are brought into main and adapted to
the main 3D model; his separate pipeline stays on his branch.
- **Lead insertion (probe bending and piercing).** `src/insertion/simulateInsertion.m` and
  `config/insertionParams.m` copied from his branch (algorithm: `docs/algorithms_and_plan.md` §4, his
  document, copied unchanged): 2D Monte Carlo per station; at a grazing contact the lead slides
  around the fascicle; at a head-on contact the fascicle is displaced if there is room and the hand
  force beats its stiffness, otherwise the lead pierces it if the force beats the puncture force.
  Forces **not calibrated** (puncture 2 mN from a tungsten needle in rabbit, Jensen 2001; hand force
  assumed): outcome frequencies are indicative, not probabilities. Added on main, off by default:
  `pierceModel = "probability"`, the tutor's rule (meeting 2026-11-23: small fascicles of our model
  very low probability, large ones high) with the angle of penetration added (user 2026-11-26): at
  the first contact with a fascicle the lead pierces it with P = 1/(1 + (r50/r)^n) · cos(γ)^m
  (`src/insertion/pierceProbability.m`; r = √(A/π), γ = angle between the heading and the perineurium
  normal), drawn once per fascicle, otherwise it bends around it; `returnToLine` steering keeps the
  exit on the line of the entry. r50 = 0.38 mm, n = 6, m = 2 `[assumed]`, chosen so that on UA01
  head-on F10 (0.17 mm) ≈ 1 %, F7 (0.28 mm) ≈ 13 %, F3 (0.43 mm) ≈ 69 %, F1 (0.69 mm) ≈ 97 %, halved at
  45°; law and values to ask the tutor (the angle dependence could be fitted to Sergi, Jensen &
  Yoshida, arXiv:1308.0525, not yet read).
  the tutor). Runs on our patient through `src/insertion/insertionGeometry.m` and
  `insertionStations.m`; study `scripts/run_insertion_3d.m` (both models, P(pierce | contact) vs
  fascicle radius); test `tests/testInsertion3D.m`. Irregular fascicle shapes
  (`src/geometry/randomFascicleShapes.m`, `nerveBoundaryRadius.m`, copied) are used in that study only.
- **Not yet in the 3D geometry.** The FEM model still has straight shanks that pierce what they
  cross (D-09) and circular fascicles; D-09 and D-11 are therefore **not** revised yet. Putting the
  bent lead and the displaced/irregular fascicles into `buildNerveGeometry3D` (PDE Toolbox, layered
  extrusion) is the next step and a design choice of its own (extend our builder, or move to a Gmsh
  mesher like Marco's).
- **Mesh convergence.** `main_pipeline_3d` split into `src/main/build3DModel.m` (patient → assembly)
  and `src/main/solve3DContacts.m` (solve, per-contact cache), with the same computation (the cache
  key checks it). `scripts/run_convergence_3d.m` follows Marco's algorithm (all element sizes scaled
  by one factor, ratio 1.3; observed order, Richardson, GCI with order check) plus the checks of
  `design.md` §4: a remeshed production level as noise floor and quantities at fixed points (R, V at
  fascicle centres, activating function).
- **Edge refinement.** Marco found R rising by 2–5 % with 5 µm elements on the pad edges (converged
  below ~5 µm). Ported as `p.mesh3d.hEdge_m` (PDE Toolbox `Hedge` on the contact edges), **0 = off**
  in production until the convergence study measures it on our mesh (turning it on changes the mesh,
  so all contacts must be solved again). pcg now accepts a stagnation below relres 1e-8 (his fix).
Not taken (on his branch, to discuss with him): the Gmsh personal pipeline (`mesh_personal.py`,
`personalDevice3D`, quicklook files), PDMS lead body and sleeve (needs the geometry step above), the
patient/function layer (`src/patient`, motor patches, V001.json; useful later for roadmap item 6),
his viewer and report changes, and his `main_pipeline_3d.m`, which is based on an older main (it
would undo D-10, D-11 and the per-contact cache).
Why / source: user 2026-11-26 ("mesh convergence algorithm, the way the probe curves in the nerve,
the probability of piercing the fascicles"); tutor meeting 2026-11-23.

## D-14 — Gmsh model (Marco's) with our solver, mesh cache, edge refinement, pcg, all 48 contacts  (2026-11-26, phases 1, 3, 4)
Decision (user 2026-11-26, after the review of Marco's push 274b485):
- **Mesher: Marco's Gmsh model** (`p.mesh3d.mesher = "gmsh"`; `src/mesh/meshPersonalDevice.m` +
  `mesh_personal.py`, copied unchanged with their dependencies: `personalDevice3D`, `quicklook3DParams`,
  `insertionScenario`, `cuffContactPlane`, `cuffDeformedRadius`, `checkCuffMeshContact`, `src/patient/`).
  The geometry gains what the tutor asked for: irregular fascicles; one random lead insertion per
  station (PDMS lead body 0.10 × 0.10 mm bent around or through the fascicles, pads placed along it,
  displaced fascicles fading out over 0.4 mm along z); an oval PDMS sleeve (0.15 mm, insulating)
  instead of the void cuff wall; flat cuff pads on the oval with a geometric indentation (3 µm extra
  closure; checked: each pad shares exactly 0.040 mm² with the epineurium). Pierce rule of the
  production insertions: the tutor's probability by fascicle size and angle (D-13,
  `pierceModel = "probability"`). **Our solver is kept**: `build3DModel` assembles the same P1 system
  with the perineurium thin layer (interface = nodes shared by endoneurium and epineurium; grounded
  nodes found by coordinates; PDMS conductivity `p.electrode.pdmsSigma_Spm` = 1e-7 S/m); then
  `solve3DContacts`, the per-contact cache, the block study and the viewer export follow as before.
  The previous PDE Toolbox model stays available (`mesher = "pdetoolbox"`). For the Gmsh model this
  revises D-03 (layered extrusion), D-05 (void cuff wall → PDMS sleeve), D-07 (no shank body → PDMS
  lead body), D-09 ("always pierce" → probability) and D-11 (rigid pads pressed in → Marco's contact plane).
- **Mesh cache** (`p.mesh3d.cacheMode`): Gmsh is not deterministic (multithreaded), and a new mesh
  invalidates every solved contact. The mesh is saved in `results/mesh/<id>_<mesher>_<key>.mat`, key =
  MD5 of anatomy, insertions, contact poses, domain and element sizes. "reuse" (default) loads it if
  the key matches, otherwise generates and saves; "regenerate" always meshes again (debugging the
  generator), keeping the previous file as `*_previous.mat` for comparison. Changing the patient or
  the fascicle layout (toy experiments) changes the key, so a new mesh is made automatically. The
  solved contacts are cached per mesh (`results/main3d/<id>_contacts/<mesh file>/`).
- **Patients** (`p.patient.source`): "generator" (our makePatient, D-09; outlines irregular with the same
  centres and areas, or circles: `p.patient.fascicleShape`), "virtual" (Marco's V<seed>), "file" (JSON,
  `docs/patient_data.md`). Several patients and toy experiments = changing id/seed/source; outputs,
  meshes and caches are all named after the patient id (user: run everything for different patients
  and toy experiments, with useful statistics as output).
- **Edge refinement on**: `p.mesh3d.hEdge_m` = 5 µm (Marco: R +2 % internal, +5 % cuff pads, converged
  below 5 µm). **pcg** solver (`p.solve3d.solver`), with automatic fallback to the direct solver if it
  does not converge (relres > 1e-8); the current balance is still checked per contact.
- **All 48 contacts solved** (`p.solve3d.sections = []`): each station has its own insertion, so
  station +k is no longer the mirror of −k and the mirror trick of D-12 does not apply (the export
  mirrors only unsolved contacts, so it simply does nothing). The block study's fibre lines follow the
  irregular outlines (`sampleFibreProfiles`).
Not done yet: the viewer still draws circular fascicles and no lead body (`export_viewer_data` only
patched so that it runs); Marco's viewer has polygon shells for this. Domain kept at 70 mm (D-02);
Marco used 60 mm.
Consequences: everything must be solved again (new mesh). Requires Python with `gmsh`
(`python -m pip install gmsh`, or `NERVE_GMSH_PYTHON` pointing to a python.exe that has it). One
insertion draw per station is one random sample: statistics need several seeds and patients.

## D-15 — Configuration selection and block: AF search over the 48 unit solutions, MRG verification, fascicle-level outcomes  (2026-11-26, phases 5, 6, block study)
Decision (user, 2026-11-26; full plan and order in `design.md` §5.6):
- **Search with the activating function, verify with MRG.** Any configuration is a weighted sum of the
  48 unit solutions (linearity, D-07), so the AF along every fibre is linear in the contact currents.
  Configurations are searched with the AF (a linear program over the signed currents, Python/scipy);
  every candidate is verified with the MRG fibres in NEURON (`src/block/`, D-12).
- **Outcomes per fascicle.** The device cannot stimulate or record single fibres, only whole fascicles
  (group, 2026-11-26). Model fibres are a population per fascicle; results, targets and constraints are
  fascicle recruitment (or block) fractions, with the 10%-of-fibres criterion (Zelechowski 2020).
  A function mixed with others inside a fascicle cannot be targeted alone.
- **Monopolar MRG threshold table first** (contact × fascicle, µA and nC): it replaces heatmap B
  (`design.md` §3) and calibrates the AF threshold on our own fibre model, instead of an external
  source (Peterson 2011, not in `Papers/`).
- **Block: one-directional distal stimulation (collision) is the next scheme to test**, as task A with
  extra constraints (same fascicle as the proximal site, distal anode stopping the downward branch),
  after a re-run of the D-12 first round on the Gmsh solution.
Alternatives: optimise directly with MRG (exact, but every evaluation is a NEURON run: too slow for a
search over 48 currents); rank contacts with heatmap B only (no thresholds, no multi-contact
patterns); kilohertz block (needs the frequency-dependent tissue model, `design.md` §5.5: kept as a
discussed alternative).
Why / source: linearity of the volume conductor (D-07); AF as excitation proxy and its limits
(course slides 02; the threshold also depends on field shape and pulse width, hence the MRG check);
first block round (`block.md` §6): collision with the same pattern and single-contact block fail.
Still open (to add here when decided): fibre diameters (one 10 µm or a distribution, source needed);
selectivity-index definition (source needed); optimiser limits (maximum active contacts, local
return or ground, contact material for the charge limit); block goal (only the fibres fired
proximally, or the whole nerve — tutor).
Consequence / limitation: AF-optimal is not MRG-optimal (a few verification loops); results depend on
the fibre population chosen inside each fascicle (positions, diameters); one insertion per station is
one random sample (robustness over insertions and patients = step A4).
