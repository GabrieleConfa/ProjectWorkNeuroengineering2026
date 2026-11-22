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
