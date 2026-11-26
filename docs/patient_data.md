# Patients: virtual or patient-specific (how to add real data)

The model runs on a **patient**: the nerve section, the outline of every fascicle, what each fascicle carries (function)
and where its motor fibres run (motor patches). A virtual patient and a measured patient have the same structure, so the
mesh, the solver, the figures and the viewer do not know the difference. Everything optional that is missing is marked
`unknown` and the model keeps working.

## 1. Choosing the patient (config/quicklook3DParams.m)

```matlab
p.patient.name = "";                   % "S01" loads data/patients/personal/S01.json (a specific patient: just drop the file there)
p.patient.file = "";                   % or the full path of any JSON file
p.patient.virtualSeed = 1;             % both empty: virtual patient number 1 = the reference patient V001 used so far
p.patient.functionScenario = "none";   % "none" (default): functions unknown, nothing assigned. Test-only hypotheses: "random", "forearmRule", "mrAIN"
```
The same options can be passed to a single call: `q = quicklook3DParams(struct('name',"S01"))`.

* **One patient or several:** the model runs on one patient at a time. To compare patients, change `p.patient.name` and run
  again. The result files are named after the case (`Cuff60`): copy `results/` somewhere before running the next patient.
* **Adding a patient** never requires touching the code: write the JSON file (section 3), `validatePatient` tells what to fix.
* Files in `data/patients/personal/`: `V001.json` (the reference virtual patient, identical to the generated one: used by the
  tests), `TEMPLATE_measured.json` (a placeholder with three fascicles to copy from).
* Virtual patients with seeds other than 1 exist (random layout by repulsion, `virtualPatient`) but are experimental: very
  dense layouts can make the shape generator or the lead insertion fail with an explanatory message.

## 2. File format (JSON, units: metres)

| Field | Meaning | Required |
|---|---|---|
| `id`, `source` | name; `"measured"` or `"virtual"` | id yes |
| `units` | must be `"m"` | yes |
| `nerve.outline` | closed outline of the nerve, list of `[x, y]` (or `nerve.semiAxes_m` and `nerve.rotation_deg`) | yes |
| `tPeri_m` | perineurium thickness (14e-6 used so far, canine data) | yes |
| `fascicles[k].outline` | closed outline of the endoneurium, list of `[x, y]`, at least 8 points | yes |
| `fascicles[k].function` | one group id (see below): the whole fascicle carries that group | optional |
| `fascicles[k].composition` | `{"AIN": 0.4, "DIGITAL": 0.6, ...}`: fraction of the fibres of each group (sums to 1) | optional |
| `fascicles[k].motorFraction` | fraction of motor fibres (about 0.13 on average in the literature) | optional |
| `fascicles[k].patches` | list of outlines of the motor patches | optional |
| `orientation` | which anatomical direction the x and y axes are (needed for the spatial hypotheses) | optional |

Functional groups (`config`: `src/patient/functionCatalog.m`): `PT` pronator teres, `FCRPL` flexor carpi radialis and palmaris
longus, `FDS` flexor digitorum superficialis, `AIN` anterior interosseous (FDP, FPL, pronator quadratus), `THENAR` thenar and
lumbrical muscles, `PCB` palmar cutaneous branch (sensory), `DIGITAL` digital sensory fibres.

## 3. From a segmented section to a patient

```matlab
nerve  = [...];                      % Nx2, metres (exported from QuPath, ImageJ, ...)
fasc   = {[...], [...], ...};        % one Nx2 matrix per fascicle
funcs  = ["", "AIN", "", ...];       % function of each fascicle ("" = unknown)
pt = patientFromOutlines("S01", nerve, fasc, 14e-6, funcs);
validatePatient(pt, quicklook3DParams());   % explains what to fix
savePatient(pt, "data/patients/personal/S01.json");
```
Then set `p.patient.file` and run `start_nerve`. The lead insertions are drawn on this patient (the section must leave room
for a lead that touches the wall on both sides: otherwise the run stops with a message that says so).

## 4. What is data and what is hypothesis

* Outlines, areas: measured for a measured patient; random (with the aggregate statistics of Pusnik 2025) for a virtual patient.
* Function: **no published map exists for the mid-arm section** (see `function_map_evidence.md`). For a virtual patient the
  composition is a stated hypothesis; `functionSource` in the patient says `assumed`, `measured` or `unknown`.
* Motor patches do not change the electrical conductivity (the endoneurial conductivity is an effective value that already
  includes the axons). They say where the motor fibres run and will be used by the recruitment calculation.
