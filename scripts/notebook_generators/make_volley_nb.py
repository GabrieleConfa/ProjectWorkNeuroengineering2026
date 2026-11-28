"""Generates notebooks/block_volley.ipynb (run once; the notebook is the deliverable)."""
import json
import sys
from pathlib import Path

OUT = Path(sys.argv[1])
cells = []


def md(s):
    """Add a markdown cell (the explanation the notebook reader sees)."""
    cells.append(("markdown", s.strip("\n")))


def code(s):
    """Add a code cell (run by the user in Jupyter)."""
    cells.append(("code", s.strip("\n")))


# =====================================================================================================
md(r'''
# The travelling volley and its anodal block

**Question (tutor, 2026-11-27).** When the proximal site fires a fascicle, its action potentials (APs)
travel along the nerve in both directions. The ones going towards the hand (the *antidromic volley*
for the sensory fibres, the motor volley for the motor ones) reach the reinnervated muscle: a twitch
the patient feels as discomfort, and noise on the EMG that controls the prosthesis. The tutor's
scheme: the **same contacts on the distal site, driven as an anode**, create a hyperpolarised region
that stops the volley. **Some single APs may still get through; what matters is that the global
activation of the target fascicle is stopped.**

This notebook

1. finds the strongest proximal stimulus that is still selective for the target fascicle (all other
   fascicles below the 10 % recruitment criterion);
2. simulates every fibre of the nerve (MRG model in NEURON) and **shows the volley travelling** along
   the nerve: one fibre (membrane voltage against position and time) and each fascicle as a whole
   (fraction of its fibres carrying an AP at each position, animated);
3. adds the **distal anode on the mirror contacts** and sweeps its current and pulse shape, first on
   the target fascicle, then on the whole nerve for the best settings;
4. judges each setting at **fascicle level** with the criterion of `docs/decisions.md` D-16 (below),
   and checks the charge against the safety limits;
5. animates the volley with the block: it travels from the proximal site and dies at the anode.

**Before running.** The fibre profiles must come from the current 3D solution (Gmsh model, all 48
contacts): run `matlab -batch "setup_paths; export_block_profiles"` from the project root first (step
B0 of `docs/design.md` §5.6). The notebook stops with a message if the profiles are older than the
solution. Then *Run All*: about 20–30 min on 7 processes (the slow cells print their time).

**Outputs:** `results/block/volley/` (figures, GIF animations, CSV tables, `summary.json`).
''')

md(r'''
## 1. Background

### 1.1 What happens after the firing: the moving AP and the volley

A fibre excited under the proximal contact fires at the node of Ranvier closest to it. From there the
AP **propagates in both directions** at the conduction velocity (≈ 50 m/s for a 10 µm MRG fibre, i.e.
≈ 1 mm every 20 µs): the depolarised node drives current through the axoplasm into its neighbours,
brings them to threshold, and so on (saltatory conduction). Behind the AP the membrane is refractory
for about a millisecond, which is why an AP cannot turn back.

A fascicle is a population of fibres. The fibres closest to the contact fire first and at the lowest
current; at the "maximum selective" current used here only part of the target fascicle fires. All the
fired fibres have the same diameter in this model, hence the same velocity: their APs travel together
as a **volley**, a compact bump of activity moving along the nerve. Its height (fraction of the
fascicle's fibres active at that position) is what a recording electrode around the fascicle would
sum up, and what the muscle "receives" when the volley arrives.

### 1.2 Anodal block, and why some APs pass

An anode raises the extracellular potential around it. Along a fibre that passes by, the activating
function (second spatial difference of the extracellular potential) is **negative under the anode**
(the membrane is hyperpolarised) and **positive on both sides** of it (the "virtual cathodes", where it
is depolarised). When the volley arrives at a strongly hyperpolarised node, the axial current from the
AP is no longer enough to bring that node to threshold, and the AP dies there: **anodal block**.

Three things limit it (all of them show up below):

* **distance**: the hyperpolarisation decays with the distance from the contact, so the fibres far
  from the anode need more current to be blocked than the near ones. The block is therefore graded
  across the fascicle: at a given current, the near fibres are blocked, the far ones are not. This is
  the "some APs pass" of the tutor;
* **virtual cathodes**: the same current that blocks near the anode depolarises the fibre on both
  sides of it. If it fires the fibre there ("anodal excitation"), the anode itself sends an AP to the
  muscle (and one to the brain). In the second block round (`docs/block.md` §5) a single anodal contact
  fired the closest fibre already at ~2 × the stimulation current, before blocking it. A slowly
  **rising** anodal current gives the membrane time to adapt (accommodation) instead of firing;
* **anode break**: when the anode is switched off abruptly, the hyperpolarised membrane rebounds and
  can fire. A slowly **falling** current avoids it. Hence the trapezoidal pulses swept below
  ("quasi-trapezoidal" pulses: van den Honert & Mortimer 1979; Ungar et al. 1986).

The anode must be on while the volley passes under it: its plateau must contain the arrival time
(distance between the sites / conduction velocity ≈ 0.4 ms after the proximal pulse).

### 1.3 When is the volley "blocked"? (D-16)

The literature on anodal block does not ask for 100 %: Tovbis et al. 2026 (FEM + MRG, the same
approach as ours, validated in the rat) measure the **anodic block = fraction of the activated fibres
whose distal-going AP is stopped**, and their best human-size cuff designs reach ≈ 90 %, never 100 %.
The tutor said the same thing. What matters for the patient is the effect on the muscle and on the
sensation, so the criterion here is at **fascicle level**, with the same 10 % yardstick that defines an
activated fascicle in the whole project (Zelechowski 2020: threshold = 10 % of the fibres recruited):

* **C0 – the sensation is kept** (added 2026-11-27): the target fascicle still delivers its APs to the
  brain end, at least 90 % of the fibres it does with the proximal pulse alone (a distal pattern could also
  hyperpolarise the brain side of the fibres);
* **C1 – the target fascicle's volley is stopped**: fewer than 10 % of the fibres of the target
  fascicle deliver an AP to the muscle end. By the same rule that says "the fascicle is activated" at
  10 %, the residual volley that reaches the muscle is below that activation level;
* **C2 – no other fascicle reaches the muscle**: in every other fascicle, fewer than 10 % of the fibres
  deliver an AP to the muscle end (the anode must not create its own volleys there);
* **C3 – no new sensation**: in every other fascicle, fewer than 10 % of the fibres send *extra* APs to
  the brain end compared with the proximal pulse alone (APs fired by the anode's virtual cathodes
  travel to the brain too);
* reported alongside: the **block ratio** of Tovbis et al. (blocked / fired, target fascicle), the
  charge per phase and the charge density on the anode, against the Shannon criterion (k = 1.75, the
  value Tovbis et al. use for platinum) and the charge-injection capacity of Pt (course slides 08; the
  CIC of gold, the material of the provided electrode, is still missing).

**Assumption to discuss with the tutor:** the 10 % level was validated on *sensory* thresholds; using it
for the residual volley at the muscle means assuming that a twitch below the fascicle's threshold
recruitment is not felt. It is a parameter (`BLOCK_LIMIT`): change it and re-run the last sections.
''')

md(r'''
## 2. Setup

Imports, project folders, NEURON, and the project's own modules:

* `mrg_fibre` – the MRG fibre in NEURON (the extracellular potential of the FEM is imposed on every
  compartment);
* `fibre_jobs` – the batch runner: one fibre simulation = one "job", sent to several processes;
  `ap_times_job` returns the AP times at every node, which is what the volley pictures need;
* `volley` – turns the per-fibre AP times into fascicle-level quantities (volley maps, reach profile,
  per-fascicle outcomes, Shannon check, GIF animations).
''')

code(r'''
# --- Standard library and the scientific stack -----------------------------------------------------
import json
import math
import os
import sys
import time
from pathlib import Path

import h5py            # FEM potential profiles written by MATLAB
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import scipy.io        # data/electrode/sites.mat
from IPython.display import Image, display

# --- Project root (Jupyter starts in notebooks/) ----------------------------------------------------
ROOT = Path.cwd().resolve()
if not (ROOT / "src" / "block").exists():
    ROOT = ROOT.parent
assert (ROOT / "src" / "block" / "mrg_fibre.py").exists(), \
    f"Cannot find the project root from {Path.cwd()}: open the notebook from the notebooks/ folder."

# --- NEURON -----------------------------------------------------------------------------------------
try:
    from neuron import h
except ImportError as err:
    raise ImportError("NEURON is not importable from this kernel. On Windows install it with the official "
                      "installer (docs/block.md §4), then start Jupyter from a NEW terminal.") from err
print("NEURON:", h.nrnversion())

# --- Project modules (src/block/) -------------------------------------------------------------------
sys.path.insert(0, str(ROOT / "src" / "block"))
from mrg_fibre import MRG_TABLE
import fibre_jobs
import volley

MECH_DIR = ROOT / "results" / "mrg" / "model"      # compiled AXNODE mechanism
BLOCK_DIR = ROOT / "results" / "block"              # inputs written by scripts/export_block_profiles.m
VOLLEY_DIR = BLOCK_DIR / "volley"                   # outputs of this notebook: shared files (threshold
VOLLEY_DIR.mkdir(parents=True, exist_ok=True)       # cache) here, each anode configuration in a subfolder
print("project root:", ROOT)
print("outputs     :", VOLLEY_DIR, "(one subfolder per anode configuration, set in section 4)")
''')

md(r'''
## 3. Parameters

Shared values come from `config/block_params.m` (through the JSON copy written by the export). The
knobs of this notebook are below; every one has a comment. Units: SI in the code, µA / ms / mm only in
labels and tables.

* **Proximal pattern**: `{contact: weight}`, cathodic < 0. AS 2 is the contact of the first two block
  rounds. `TARGET_FASCICLE = None` lets the notebook pick the fascicle the pattern excites most easily
  (lowest 10 % threshold): with the Gmsh model every station has its own random lead insertion, so the
  fascicle around each contact may differ from the old model.
* **Maximum selective current**: `I_prox = SEL_MARGIN × (lowest 10 % threshold among the other
  fascicles)`: just below the current at which a second fascicle would count as activated.
* **Anode**: plateau from `PLATEAU_START_s` to `PLATEAU_END_s` after the proximal onset (it must contain
  the arrival of the volley, checked in section 5); ramps of `RISE_TIMES_s` (0 = rectangular); currents
  `ANODE_AMPS_A`.
* **Where the anodal current returns** (`ANODE_RETURN`, D-18). `None` = **monopolar** anode: the current
  goes to the far ground (the outer muscle surface). In the first runs this sent the whole anodal current
  along the nerve inside the insulating sleeve, and at the sleeve end (z = ±25 mm) it fired every fascicle
  before any block (`docs/block.md` §6). `"distal_site"` (default) = **balanced** anode: the same current is
  taken back by all the other contacts of the distal site (23 contacts, split equally, sum of the weights
  = 0), so almost no net current flows to the sleeve ends, and each return contact carries only 1/23 of the
  current. The return contacts act as weak cathodes during the anodal pulse: the AF check below and the
  whole-nerve check (section 7) tell whether they fire other fibres. `"station_cuff"` returns through the 2
  cuff pads of the anode's station only (the static AF check predicts that they would push the nearest
  fascicle much harder, D-18); a dict `{contact: weight}` with negative weights chooses other contacts.
* **Criteria** (section 1.3): `BLOCK_LIMIT` for C1, `SIDE_LIMIT` for C2 and C3.
''')

code(r'''
PATIENT = "UA01"
P = json.loads((BLOCK_DIR / f"{PATIENT}_block_params.json").read_text())
MRG = P["mrg"]                         # MRG settings: diameter, length, time step, AP detection...

# --- Proximal stimulation ---------------------------------------------------------------------------
PROX_PATTERN = {2: -1.0}               # {contact id: weight}, cathodic < 0
TARGET_FASCICLE = None                 # None = automatic (lowest 10 % threshold for PROX_PATTERN)
FIBRE_DIAM_um = MRG["fibreDiam_um"]    # 10 um (config); allowed: keys of MRG_TABLE
RECRUIT = P["recruitFraction"]         # 0.10: the 10 %-of-fibres criterion (Zelechowski 2020)
SEL_MARGIN = 0.95                      # I_prox = 0.95 x the lowest 10 % threshold of the other fascicles

# --- Timing (every anode time is counted from the ONSET of the proximal pulse) -----------------------
T_TH_ON_s    = P["pulse"]["stimOnset_s"]   # pulse onset in the threshold searches (0.5 ms of rest)
T_PROX_ON_s  = 1.5e-3                  # proximal onset in the block runs: > the longest ramp, so that a
                                       # ramp can start BEFORE the proximal pulse with the fibre at rest
PROX_WIDTH_s = P["pulse"]["stimWidth_s"]   # 50 us (Zelechowski 2020)
T_AFTER_s    = P["pulse"]["tAfter_s"]      # simulated time after the last pulse

# --- Distal anode -----------------------------------------------------------------------------------
PLATEAU_START_s  = 0.1e-3              # full amplitude from 0.1 ms ...
PLATEAU_END_s    = 0.9e-3              # ... to 0.9 ms after the proximal onset
RISE_TIMES_s     = [0.0, 0.3e-3, 1.0e-3]   # ramp lengths (0 = rectangular)
FALL_EQUALS_RISE = True                # symmetric trapezoid (slow fall against anode break)
RAMP_STEP_s      = 10e-6               # a ramp = staircase of 10 us steps (2 MRG time steps)
ANODE_AMPS_A     = np.geomspace(10e-6, 2e-3, 12)   # 10 uA ... 2 mA, x1.6 per step
ANODE_RETURN     = "distal_site"       # None = monopolar (far ground); "distal_site" = balanced, returned by
                                       # all the other contacts of the distal site, split equally (D-18);
                                       # "station_cuff" = by the 2 cuff pads of the anode's station;
                                       # or {contact: weight < 0}
DIST_PATTERN     = None                # None: the mirror anode + ANODE_RETURN above. Otherwise the whole
                                       # distal pattern {contact: weight} (anode > 0, cathode < 0, largest
                                       # |weight| = 1), or "design:<name>" = a candidate written by
                                       # notebooks/block_anode_design.ipynb (D-19); ANODE_RETURN is then ignored
SLEEVE_HALF_m    = 25e-3               # half-length of the insulating sleeve (p.cuff.length_m / 2 in
                                       # config/pw6_3d_params.m): its ends are marked in the AF check

# --- Fascicle-level criteria (D-16) -----------------------------------------------------------------
BLOCK_LIMIT = RECRUIT                  # C1: target fibres delivering an AP to the muscle < 10 % of the fascicle
SIDE_LIMIT  = RECRUIT                  # C2/C3: other fascicles < 10 % at the muscle end / extra at the brain end
N_CHECK     = 3                        # settings of section 6 re-checked on the whole nerve in section 7

# --- Safety (rough) ---------------------------------------------------------------------------------
SHANNON_K = 1.75                       # Shannon 1992 k used for Pt by Tovbis et al. 2026
CIC_PT_uC_cm2 = (20.0, 300.0)          # Pt charge-injection capacity, course slides 08 (0.02-0.3 mC/cm^2)

# --- Pictures ---------------------------------------------------------------------------------------
VOLLEY_WINDOW_s = 0.1e-3               # an AP is drawn at a node for 0.1 ms after it crosses 0 mV there
FRAME_STEP_s    = 20e-6                # animation frames every 20 us (~1 mm of travel at 50 m/s)
GIF_FPS         = 15

N_WORKERS = max(1, (os.cpu_count() or 2) - 1)

assert T_PROX_ON_s + PLATEAU_START_s - max(RISE_TIMES_s) >= 0.2e-3, \
    "the longest ramp would start before t = 0.2 ms (fibre not at rest): increase T_PROX_ON_s"
assert MRG["nodeAlignment"] == "proximalStation", \
    "the shared fibre object needs nodeAlignment = proximalStation (config/block_params.m)"
print(f"fibre D = {FIBRE_DIAM_um} um (internode {MRG_TABLE[FIBRE_DIAM_um][5]} um), {N_WORKERS} worker processes")
''')

md(r'''
## 4. Data: FEM potential along each fibre, electrode layout, anode pattern

The profile file holds, for every solved contact and every fibre line, the FEM potential along the
fibre for 1 A injected by that contact alone. Any pattern is a weighted sum of these profiles
(superposition, the FEM problem is linear).

**Anode pattern = the same electrode on the distal site** (supervisor's rule, D-12; confirmed by the
user 2026-11-27): each proximal contact is paired with the contact of the same type, radius and angle
on the opposite station (section −k → +k), and its weight changes sign (cathode → anode).
''')

code(r'''
# Load the FEM potential along every fibre line (HDF5 written by scripts/export_block_profiles.m), the
# electrode layout, and build the distal anode pattern (ANODE_RETURN / DIST_PATTERN knobs of section 3).
# Stops if the profiles are older than the 3D solution (they must be re-exported after every re-solve).
H5_PATH = BLOCK_DIR / f"{PATIENT}_fibre_profiles.h5"
FEM_MAT = ROOT / "results" / "main3d" / f"{PATIENT}_3d.mat"
if FEM_MAT.exists() and H5_PATH.stat().st_mtime < FEM_MAT.stat().st_mtime:
    raise RuntimeError(f"{H5_PATH.name} is older than {FEM_MAT.name}: the profiles come from an old solution. "
                       "Run  matlab -batch \"setup_paths; export_block_profiles\"  first (step B0).")

with h5py.File(H5_PATH, "r") as f:
    z_m         = f["z_m"][:].ravel()                        # sample positions along the nerve [m]
    fibre_x_m   = f["x_m"][:].ravel()                        # fibre line positions in the section [m]
    fibre_y_m   = f["y_m"][:].ravel()
    fascicle    = f["fascicle"][:].ravel().astype(int)       # fascicle of each fibre line
    contact_id  = f["contact_id"][:].ravel().astype(int)
    contact_z_m = f["contact_z_m"][:].ravel()
    n_contacts, n_fibres, n_z = f["V_perA"].shape
    fem_source = f.attrs["source"]
fem_source = fem_source.decode() if isinstance(fem_source, bytes) else str(fem_source)
CONTACT_Z = dict(zip(contact_id.tolist(), contact_z_m.tolist()))
FASCICLES = np.unique(fascicle)
N_PER_FASC = {int(k): int(np.sum(fascicle == k)) for k in FASCICLES}
print(f"FEM solution: {fem_source}")
print(f"{n_fibres} fibres in {len(FASCICLES)} fascicles {N_PER_FASC}; {n_contacts} contacts")

sites_all = scipy.io.loadmat(str(ROOT / "data" / "electrode" / "sites.mat"))["sites"]
sites_all = pd.DataFrame(sites_all, columns=["AS", "type", "section", "z_mm", "r_mm", "theta_deg", "phi_deg",
                                             "thickness_mm", "sideZ_mm", "sideXY_mm"])
sites_all["AS"] = sites_all["AS"].astype(int)
sites_all = sites_all.set_index("AS")


def contact_area_m2(c):
    """Area of one face of contact c [m^2] (sites.mat: sideZ x sideXY), used for the charge density."""
    s = sites_all.loc[c]
    return float(s.sideZ_mm * s.sideXY_mm) * 1e-6


def layout_partners(sites_all):
    """Proximal contact -> the same electrode on the distal site (same type, radius, angle; section -k -> +k)."""
    out = {}
    for c, s in sites_all[sites_all.section < 0].iterrows():
        cand = sites_all[(sites_all.section == -s.section) & (sites_all.type == s.type) &
                         (np.isclose(sites_all.theta_deg, s.theta_deg)) &
                         ((sites_all.r_mm.isna() & np.isnan(s.r_mm)) | np.isclose(sites_all.r_mm, s.r_mm))]
        if len(cand) == 1:
            out[int(c)] = int(cand.index[0])
    return out


MIRROR = layout_partners(sites_all)
missing = [c for c in PROX_PATTERN if c not in MIRROR or MIRROR[c] not in CONTACT_Z]
assert not missing, f"no solved distal partner for contacts {missing} in the profile file"
ANODE_PATTERN = {MIRROR[c]: -w for c, w in PROX_PATTERN.items()}
assert all(w > 0 for w in ANODE_PATTERN.values()), "the proximal weights must all be cathodic (< 0)"
PROX_NAME = "+".join(f"AS{c}" for c in PROX_PATTERN)
ANODE_NAME = "+".join(f"AS{c}" for c in ANODE_PATTERN)
Z_PROX_m = float(np.mean([CONTACT_Z[c] for c in PROX_PATTERN]))
Z_DIST_m = float(np.mean([CONTACT_Z[c] for c in ANODE_PATTERN]))
assert Z_PROX_m < 0 < Z_DIST_m, "proximal contacts at z < 0, distal at z > 0 (D-12)"

# --- Return of the anodal current (D-18): build the full distal pattern ------------------------------
ANODE_CONTACTS = dict(ANODE_PATTERN)                       # the mirror anode(s) alone
if ANODE_RETURN is None:
    RETURN_PATTERN, ANODE_TAG = {}, f"monopolar_{ANODE_NAME}"
elif ANODE_RETURN == "distal_site":
    others = [c for c in CONTACT_Z if CONTACT_Z[c] > 0 and c not in ANODE_CONTACTS]
    RETURN_PATTERN = {c: -sum(ANODE_CONTACTS.values()) / len(others) for c in sorted(others)}
    ANODE_TAG = f"balanced_{ANODE_NAME}_ret-distal-site"
elif ANODE_RETURN == "station_cuff":
    secs = sites_all.loc[list(ANODE_CONTACTS)].section.unique()
    pads = [int(c) for c in sites_all[sites_all.section.isin(secs) & (sites_all.type != 1)].index
            if int(c) in CONTACT_Z and int(c) not in ANODE_CONTACTS]
    assert pads, f"no solved cuff pad on the anode's station {secs}"
    RETURN_PATTERN = {c: -sum(ANODE_CONTACTS.values()) / len(pads) for c in pads}
    ANODE_TAG = f"balanced_{ANODE_NAME}_ret" + "-".join(str(c) for c in pads)
else:
    RETURN_PATTERN = {int(c): float(w) for c, w in ANODE_RETURN.items()}
    assert all(w < 0 for w in RETURN_PATTERN.values()), "return weights must be negative (cathodic)"
    ANODE_TAG = f"custom_{ANODE_NAME}_ret" + "-".join(str(c) for c in RETURN_PATTERN)
ANODE_PATTERN = {**ANODE_CONTACTS, **RETURN_PATTERN}        # the distal channel actually driven
if DIST_PATTERN is not None:                               # a designed / custom pattern replaces all of it
    if isinstance(DIST_PATTERN, str) and DIST_PATTERN.startswith("design:"):
        dname = DIST_PATTERN.split(":", 1)[1]
        cands = json.loads((BLOCK_DIR / "anode_design" / "candidates.json").read_text())
        ANODE_PATTERN = {int(c): float(w) for c, w in cands[dname]["weights"].items()}
        ANODE_TAG = f"design_{dname}"
    else:
        ANODE_PATTERN = {int(c): float(w) for c, w in DIST_PATTERN.items()}
        ANODE_TAG = "custom_" + "_".join(f"{c}{'p' if w > 0 else 'm'}" for c, w in ANODE_PATTERN.items())
    ANODE_CONTACTS = {c: w for c, w in ANODE_PATTERN.items() if w > 0}
    RETURN_PATTERN = {c: w for c, w in ANODE_PATTERN.items() if w < 0}
    assert ANODE_CONTACTS, "the distal pattern needs at least one anode (weight > 0)"
    Z_DIST_m = float(np.average([CONTACT_Z[c] for c in ANODE_CONTACTS], weights=list(ANODE_CONTACTS.values())))
    ANODE_NAME = "+".join(f"AS{c}" for c in ANODE_CONTACTS) + ("/" + "+".join(f"AS{c}" for c in RETURN_PATTERN)
                                                              + "(-)" if RETURN_PATTERN else "")
elif RETURN_PATTERN:
    ANODE_NAME += ("/rest of distal site(-)" if ANODE_RETURN == "distal_site"
                   else "/" + "+".join(f"AS{c}" for c in RETURN_PATTERN) + "(-)")
OUT_DIR = VOLLEY_DIR / f"{PROX_NAME}_{ANODE_TAG}"
OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"proximal cathode {PROX_PATTERN} at z = {Z_PROX_m*1e3:.1f} mm  ->  distal pattern {ANODE_PATTERN} "
      f"at z = {Z_DIST_m*1e3:.1f} mm (net current to the far ground: {sum(ANODE_PATTERN.values()):+.2f} x I)")
print("outputs of this configuration:", OUT_DIR)
''')

md(r'''
## 5. Proximal stimulation and the volley without block

### 5.1 Fibre model and pulse helpers

One MRG fibre object is built in every process (a node of Ranvier exactly at the proximal station) and
re-used for every fibre line: the lines differ only by the extracellular potential applied to them.

`mrg_fibre` applies rectangular pulses; a trapezoid is written as a staircase of 10 µs steps.
''')

code(r'''
# Build the MRG fibre (one shared NEURON fibre, nodes aligned under the proximal contact) and the helpers
# that turn a current into stimulator drives: rect (one pulse), the proximal drive, the anode pulses with
# optional linear ramps written as a staircase of short pulses.
t0 = time.time()
runner = fibre_jobs.FibreRunner(H5_PATH, MRG, FIBRE_DIAM_um, Z_PROX_m, MECH_DIR, N_WORKERS)
FIBRE = runner.fibre
NODE_B, NODE_M = runner.node_brain, runner.node_muscle
N_NODES = FIBRE.n_nodes
print(f"MRG fibre: {N_NODES} nodes, a node every {FIBRE.internode_m*1e3:.2f} mm; APs counted at "
      f"z = {FIBRE.z_node_m[NODE_B]*1e3:.1f} mm (brain end) and z = {FIBRE.z_node_m[NODE_M]*1e3:.1f} mm "
      f"(muscle end)  [{time.time()-t0:.1f} s]")


def rect(t_on_s, width_s, amp_A):
    """One rectangular pulse as the (t_on, t_off, amplitude) list mrg_fibre expects."""
    return [(t_on_s, t_on_s + width_s, amp_A)]


def trapezoid(t_plateau_on_s, t_plateau_off_s, amp_A, rise_s, fall_s, step_s=RAMP_STEP_s):
    """Trapezoid as a staircase: the rise ENDS at t_plateau_on_s, the fall STARTS at t_plateau_off_s."""
    pulses = []
    n = int(round(rise_s / step_s))
    if n > 0:
        edges = t_plateau_on_s - rise_s + np.arange(n + 1) * rise_s / n
        pulses += [(edges[k], edges[k + 1], amp_A * (k + 0.5) / n) for k in range(n)]
    pulses.append((t_plateau_on_s, t_plateau_off_s, amp_A))
    n = int(round(fall_s / step_s))
    if n > 0:
        edges = t_plateau_off_s + np.arange(n + 1) * fall_s / n
        pulses += [(edges[k], edges[k + 1], amp_A * (1 - (k + 0.5) / n)) for k in range(n)]
    return pulses


def charge_C(pulses):
    """Charge of one phase [C] (time integral of the current)."""
    return sum((off - on) * amp for on, off, amp in pulses)


def end_s(pulses):
    return max(off for _, off, _ in pulses)


def start_s(pulses):
    return min(on for on, _, _ in pulses)


def prox_drive(amp_A, t_on_s=T_PROX_ON_s):
    """Proximal channel: the pattern with one rectangular cathodic pulse."""
    return (PROX_PATTERN, rect(t_on_s, PROX_WIDTH_s, amp_A))


def anode_pulses(amp_A, rise_s):
    """Distal anodal pulse: plateau PLATEAU_START..PLATEAU_END after the proximal onset."""
    fall_s = rise_s if FALL_EQUALS_RISE else 0.0
    return trapezoid(T_PROX_ON_s + PLATEAU_START_s, T_PROX_ON_s + PLATEAU_END_s, amp_A, rise_s, fall_s)


def t_stop_for(*pulse_lists):
    return max(end_s(p) for p in pulse_lists) + T_AFTER_s


def anode_safety(pulses):
    """Charge and density per phase on the distal contact with the highest density, Shannon, CIC."""
    worst = max(ANODE_PATTERN, key=lambda c: abs(ANODE_PATTERN[c]) / contact_area_m2(c))   # highest density
    q = charge_C(pulses) * abs(ANODE_PATTERN[worst])
    ok, Q_uC, D = volley.shannon_ok(q, contact_area_m2(worst), SHANNON_K)
    return {"charge_nC": Q_uC * 1e3, "density_uC_cm2": D, "shannon_ok": ok, "worst_contact": int(worst),
            "below_Pt_CIC_low": D <= CIC_PT_uC_cm2[0], "below_Pt_CIC_high": D <= CIC_PT_uC_cm2[1]}


T_STOP_BASE_s = T_PROX_ON_s + PROX_WIDTH_s + T_AFTER_s
''')

md(r'''
### 5.2 Thresholds of every fibre and the maximum selective current

For every fibre of the nerve, the smallest current of the proximal pattern (one 50 µs cathodic pulse)
that fires it (bisection, `fibre_jobs.threshold_job`). Sorting the thresholds of a fascicle gives its
**recruitment curve**; its 10 % threshold is the current at which 10 % of its fibres fire.

* the **target** is the fascicle with the lowest 10 % threshold (unless set by hand);
* the **maximum selective current** is `SEL_MARGIN` × the lowest 10 % threshold among the *other*
  fascicles: the strongest stimulus that still activates only the target fascicle.

This is the slowest single step (~453 fibres × ~10 runs): the result is cached in `OUT_DIR` together
with the inputs it depends on, and re-used if nothing changed.
''')

code(r'''
# Threshold of every fibre for the proximal pattern (doubling + bisection, fibre_jobs.threshold_job),
# cached on disk by its inputs; then the 10 %-of-fibres threshold of every fascicle, the target fascicle
# (lowest threshold) and the maximum selective current I_prox = SEL_MARGIN x the lowest other threshold.
TH_CSV = VOLLEY_DIR / f"{PATIENT}_{PROX_NAME}_D{FIBRE_DIAM_um:g}_thresholds.csv"
TH_KEY = VOLLEY_DIR / f"{PATIENT}_{PROX_NAME}_D{FIBRE_DIAM_um:g}_thresholds_key.json"
key = {"h5_mtime": H5_PATH.stat().st_mtime, "h5_size": H5_PATH.stat().st_size,
       "pattern": {str(c): w for c, w in PROX_PATTERN.items()}, "diam_um": FIBRE_DIAM_um,
       "width_s": PROX_WIDTH_s, "t_on_s": T_TH_ON_s, "z_align_m": Z_PROX_m,
       "search": [P["searchStart_A"], P["searchMax_A"], P["searchTol"]], "mrg": MRG}

t0 = time.time()
if TH_CSV.exists() and TH_KEY.exists() and json.loads(TH_KEY.read_text()) == json.loads(json.dumps(key)):
    th_A = pd.read_csv(TH_CSV)["threshold_A"].to_numpy()
    print(f"thresholds loaded from the cache {TH_CSV.name}")
else:
    t_stop = T_TH_ON_s + PROX_WIDTH_s + T_AFTER_s
    jobs = [(i, PROX_PATTERN, T_TH_ON_s, PROX_WIDTH_s, t_stop, P["searchStart_A"], P["searchMax_A"],
             P["searchTol"]) for i in range(n_fibres)]
    print(f"threshold search on {n_fibres} fibres, {N_WORKERS} processes ...")
    th_A = np.array(runner.map(fibre_jobs.threshold_job, jobs, chunksize=2))
    pd.DataFrame({"fibre": np.arange(n_fibres), "fascicle": fascicle, "threshold_A": th_A}).to_csv(TH_CSV, index=False)
    TH_KEY.write_text(json.dumps(key, indent=1))
print(f"[{time.time()-t0:.0f} s]  fibres that never fire up to {P['searchMax_A']*1e6:.0f} uA: {np.isnan(th_A).sum()}")


def th10_A(k):
    """10 % threshold of fascicle k [A] (NaN if fewer than 10 % of its fibres fire up to searchMax)."""
    t = np.sort(th_A[fascicle == k])
    n_need = max(1, math.ceil(RECRUIT * len(t)))
    t = t[np.isfinite(t)]
    return t[n_need - 1] if len(t) >= n_need else np.nan


TH10 = {int(k): th10_A(k) for k in FASCICLES}
TARGET = TARGET_FASCICLE if TARGET_FASCICLE is not None else int(min(TH10, key=lambda k: np.inf if np.isnan(TH10[k]) else TH10[k]))
others = [TH10[k] for k in TH10 if k != TARGET and np.isfinite(TH10[k])]
I_PROX = SEL_MARGIN * min(others) if others else P["searchMax_A"]
SELECTIVE = I_PROX > TH10[TARGET]
if not SELECTIVE:
    print(f"WARNING: {PROX_NAME} is not selective for F{TARGET} (another fascicle reaches 10 % first); "
          f"falling back to {P['stimFactor']} x the target's 10 % threshold")
    I_PROX = P["stimFactor"] * TH10[TARGET]

REC = pd.DataFrame([{"fascicle": k, "n_fibres": N_PER_FASC[k], "th10_uA": TH10[k] * 1e6,
                     "th10_nC": TH10[k] * PROX_WIDTH_s * 1e9,
                     "recruited_at_I_prox": float(np.mean(th_A[fascicle == k] <= I_PROX))}
                    for k in TH10]).set_index("fascicle")
REC.to_csv(OUT_DIR / "recruitment_table.csv")
print(f"target F{TARGET}: 10 % threshold {TH10[TARGET]*1e6:.1f} uA; next fascicle at {min(others)*1e6:.1f} uA "
      f"(window x{min(others)/TH10[TARGET]:.1f})")
print(f"I_prox = {I_PROX*1e6:.1f} uA = {I_PROX*PROX_WIDTH_s*1e9:.2f} nC -> recruits "
      f"{REC.recruited_at_I_prox[TARGET]*100:.0f} % of F{TARGET}")
REC.round(2)
''')

md(r'''
### How to read this figure: recruitment curves

One curve per fascicle: the fraction of its fibres whose threshold is below the current on the x axis
(log scale). The target fascicle is the thick line. The dotted horizontal line is the 10 % criterion,
the vertical line the maximum selective current `I_prox`: where it crosses the target curve is the
recruited fraction of the target (the size of the volley that must be blocked); all the other curves
must still be below 10 % there.
''')

code(r'''
# Recruitment curves: fraction of each fascicle fired against the proximal current (from the thresholds).
amps = np.geomspace(np.nanmin(th_A) * 0.8, P["searchMax_A"], 200)
fig, ax = plt.subplots(figsize=(8, 4.5))
cmap = plt.get_cmap("tab10")
for k in FASCICLES:
    t = th_A[fascicle == k]
    frac = [(t <= a).mean() for a in amps]
    ax.semilogx(amps * 1e6, frac, color=cmap((k - 1) % 10), lw=3 if k == TARGET else 1.2,
                label=f"F{k}" + (" (target)" if k == TARGET else ""))
ax.axhline(RECRUIT, color="k", ls=":", lw=0.8)
ax.axvline(I_PROX * 1e6, color="k", lw=0.8)
ax.text(I_PROX * 1e6, 1.02, " I_prox", fontsize=8)
ax.set(xlabel=f"current of {PROX_NAME}, 50 us cathodic pulse [uA]", ylabel="fraction of the fascicle's fibres fired",
       ylim=(0, 1.08), title=f"Recruitment curves, MRG {FIBRE_DIAM_um} um")
ax.legend(fontsize=7, ncol=2)
fig.savefig(OUT_DIR / "recruitment_curves.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
### 5.3 The volley without block: every fibre at I_prox

One run per fibre with the proximal pulse alone, recording the AP times at every node. This is the
**reference** for everything below: which fibres fire, and what reaches each end. With the proximal
pulse alone every fired fibre should send one AP to the brain end and one to the muscle end.
''')

code(r'''
# Reference run: every fibre with the proximal pulse alone. AP_BASE (AP times at every node) and
# OUT_BASE (per-fascicle fractions reaching the brain and the muscle end) are the baseline of the criteria.
t0 = time.time()
AP_BASE = runner.run_ap_times([(i, [prox_drive(I_PROX)], T_STOP_BASE_s) for i in range(n_fibres)])
OUT_BASE = volley.fascicle_outcomes(AP_BASE, fascicle, NODE_B, NODE_M)
FIRED = (volley.counts_at(AP_BASE, NODE_B) + volley.counts_at(AP_BASE, NODE_M)) > 0
TGT = np.flatnonzero(fascicle == TARGET)
assert FIRED[TGT].any(), "the proximal pulse fires no fibre of the target fascicle"
print(f"[{time.time()-t0:.0f} s]  fired fibres: {FIRED.sum()} "
      f"({FIRED[TGT].sum()} of the {len(TGT)} fibres of F{TARGET}); predicted from the thresholds: "
      f"{int(np.sum(th_A <= I_PROX))}")
OUT_BASE.to_csv(OUT_DIR / "outcomes_proximal_alone.csv")
OUT_BASE[["n_fibres", "brain", "muscle", "frac_brain", "frac_muscle"]].round(3)
''')

md(r'''
**Conduction velocity and arrival time at the distal site.** From the fired fibres of the target
fascicle: the velocity between two nodes on the distal side, and the time at which the volley reaches
the node under the anode. The anode plateau must contain all the arrival times.
''')

code(r'''
# Conduction velocity of the target fibres and arrival time of the volley at the distal site: the anode
# plateau (PLATEAU_START_s ... PLATEAU_END_s) must contain every arrival.
n_a, n_b = FIBRE.node_index(Z_DIST_m - 5e-3), FIBRE.node_index(Z_DIST_m + 10e-3)
n_d = FIBRE.node_index(Z_DIST_m)
fired_tgt = TGT[FIRED[TGT]]
cv = [(FIBRE.z_node_m[n_b] - FIBRE.z_node_m[n_a]) /
      (volley.first_crossing_s(AP_BASE[i], n_b) - volley.first_crossing_s(AP_BASE[i], n_a)) for i in fired_tgt]
arr = np.array([volley.first_crossing_s(AP_BASE[i], n_d) for i in fired_tgt]) - T_PROX_ON_s
CV_mps = float(np.nanmedian(cv))
ARRIVAL_s = float(np.nanmedian(arr))
print(f"conduction velocity {CV_mps:.1f} m/s; the volley reaches the distal site "
      f"{np.nanmin(arr)*1e3:.2f}-{np.nanmax(arr)*1e3:.2f} ms after the proximal onset "
      f"({(Z_DIST_m - Z_PROX_m)*1e3:.1f} mm apart)")
if not (PLATEAU_START_s < np.nanmin(arr) and np.nanmax(arr) < PLATEAU_END_s):
    print("WARNING: the anode plateau does not contain every arrival: change PLATEAU_START_s / PLATEAU_END_s")

# frames of the animations: from just before the pulse until the volley has left the fibre
travel_s = max(abs(FIBRE.z_node_m[-1] - Z_PROX_m), abs(Z_PROX_m - FIBRE.z_node_m[0])) / CV_mps
FRAMES_BASE = np.arange(T_PROX_ON_s - 0.1e-3, T_PROX_ON_s + travel_s + 0.3e-3, FRAME_STEP_s)
SITE_MARKS = [(Z_PROX_m, f"proximal {PROX_NAME}"), (Z_DIST_m, f"distal {ANODE_NAME}")]
''')

md(r'''
### 5.4 One fibre: the moving AP

The most excitable fired fibre of the target fascicle, with the membrane voltage recorded at every
node.

**How to read this figure.** *Left*: membrane voltage (colour) at each node (y, position along the
nerve) against time (x). The AP is the bright diagonal band: it starts at the proximal site and moves
in both directions; the slope is the conduction velocity. *Right*: snapshots of the voltage along the
fibre at a few instants, the same wave seen as a profile moving along the nerve (the "moving AP"). The
GIF below shows it in motion.
''')

code(r'''
# One fibre of the target fascicle: membrane voltage at all nodes vs time (raster) and snapshots of the
# travelling AP.
I_EX = int(fired_tgt[np.nanargmin(th_A[fired_tgt])])
_, _, rec = runner.record(I_EX, [prox_drive(I_PROX)], T_STOP_BASE_s)
t_rec, v_rec = rec["t_s"], rec["v_V"]
zn_mm = FIBRE.z_node_m * 1e3

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.1, 1]})
im = a1.pcolormesh((t_rec - T_PROX_ON_s) * 1e3, zn_mm, v_rec * 1e3, cmap="magma", vmin=-100, vmax=40, shading="auto")
for zz, lab in SITE_MARKS:
    a1.axhline(zz * 1e3, color="c", lw=0.8, ls="--")
    a1.text(0.02, zz * 1e3, lab, color="c", fontsize=7, va="bottom", transform=a1.get_yaxis_transform())
a1.set(xlim=(-0.1, travel_s * 1e3 + 0.3), xlabel="time from the proximal onset [ms]",
       ylabel="node position z [mm]  (brain bottom, muscle top)",
       title=f"Fibre {I_EX} (F{fascicle[I_EX]}): membrane voltage at the nodes")
fig.colorbar(im, ax=a1, label="V_m [mV]")
for dt_ms in [0.05, 0.15, 0.3, 0.45, 0.6]:
    k = np.searchsorted(t_rec, T_PROX_ON_s + dt_ms * 1e-3)
    a2.plot(zn_mm, v_rec[:, k] * 1e3, label=f"t = {dt_ms:.2f} ms")
for zz, _ in SITE_MARKS:
    a2.axvline(zz * 1e3, color="grey", ls="--", lw=0.8)
a2.set(xlabel="z along the nerve [mm]  (brain <- | -> muscle)", ylabel="V_m at the nodes [mV]",
       title="Snapshots: two APs leave the proximal site in opposite directions")
a2.legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT_DIR / "single_fibre_moving_AP.png", dpi=130, bbox_inches="tight")
plt.show()

idx = np.clip(np.searchsorted(t_rec, FRAMES_BASE), 0, len(t_rec) - 1)
gif = volley.animate_profiles(OUT_DIR / "single_fibre_moving_AP.gif", {f"fibre {I_EX}": (v_rec[:, idx] * 1e3).T},
                              FRAMES_BASE, FIBRE.z_node_m, ylabel="V_m at the nodes [mV]", ylim=(-110, 50),
                              markers=SITE_MARKS, t_ref_s=T_PROX_ON_s, fps=GIF_FPS,
                              title=f"Fibre {I_EX} (F{fascicle[I_EX]}), proximal pulse {I_PROX*1e6:.0f} uA, no block")
display(Image(filename=str(gif)))
''')

md(r'''
### 5.5 The whole fascicle: the volley

Now the population. For each fascicle, at each instant and each node, the fraction of **all** its
fibres that carry an AP there (`volley.volley_map`; the denominator is the whole fascicle, so the
values compare directly with the 10 % criterion).

**How to read this figure.** *Left*: the volley of the target fascicle in space (y) and time (x): the
same diagonal band as for one fibre, now with the colour giving the fraction of the fascicle active
there. Its maximum equals the recruited fraction. *Right*: the volley of each fascicle as a profile
along the nerve at one instant (when the distal-going volley is half-way to the distal site). The GIF
shows the profiles moving: one bump per fascicle and per direction; the other fascicles should stay
flat (below 10 %) at the maximum selective current.
''')

code(r'''
# The whole target fascicle: fraction of its fibres carrying an AP at each node and time (volley map),
# and the volley profile of every fascicle half-way to the distal site.
VMAP_BASE = {int(k): volley.volley_map(AP_BASE, np.flatnonzero(fascicle == k), N_NODES, FRAMES_BASE,
                                       VOLLEY_WINDOW_s, n_norm=N_PER_FASC[int(k)]) for k in FASCICLES}
active = [k for k in VMAP_BASE if VMAP_BASE[k].max() > 0]
colours = {f"F{k}": cmap((k - 1) % 10) for k in FASCICLES}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw={"width_ratios": [1.1, 1]})
im = a1.pcolormesh((FRAMES_BASE - T_PROX_ON_s) * 1e3, zn_mm, VMAP_BASE[TARGET].T, cmap="viridis",
                   vmin=0, vmax=max(0.05, VMAP_BASE[TARGET].max()), shading="auto")
for zz, lab in SITE_MARKS:
    a1.axhline(zz * 1e3, color="w", lw=0.8, ls="--")
a1.set(xlabel="time from the proximal onset [ms]", ylabel="node position z [mm]",
       title=f"Volley of F{TARGET}: fraction of its fibres carrying an AP")
fig.colorbar(im, ax=a1, label="fraction of the fascicle")
k_half = np.searchsorted(FRAMES_BASE, T_PROX_ON_s + ARRIVAL_s / 2)
for k in active:
    a2.plot(zn_mm, VMAP_BASE[k][k_half], color=colours[f"F{k}"], lw=2.5 if k == TARGET else 1, label=f"F{k}")
a2.axhline(RECRUIT, color="k", ls=":", lw=0.8)
for zz, _ in SITE_MARKS:
    a2.axvline(zz * 1e3, color="grey", ls="--", lw=0.8)
a2.set(xlabel="z [mm]", ylabel="fraction of the fascicle with an AP",
       title=f"t = {(FRAMES_BASE[k_half]-T_PROX_ON_s)*1e3:.2f} ms: volleys along the nerve")
a2.legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT_DIR / "volley_no_block.png", dpi=130, bbox_inches="tight")
plt.show()

gif = volley.animate_profiles(OUT_DIR / "volley_no_block.gif", {f"F{k}": VMAP_BASE[k] for k in active},
                              FRAMES_BASE, FIBRE.z_node_m, ylabel="fraction of the fascicle with an AP",
                              ylim=(0, max(0.15, max(VMAP_BASE[k].max() for k in active)) * 1.1),
                              markers=SITE_MARKS, t_ref_s=T_PROX_ON_s, fps=GIF_FPS, colours=colours,
                              title=f"Volleys after the proximal pulse ({I_PROX*1e6:.0f} uA on {PROX_NAME}), no block")
display(Image(filename=str(gif)))
''')

md(r'''
## 6. The distal anode on the target fascicle: current × pulse shape

For every setting (ramp length × anodal current) and every fibre of the target fascicle, two runs:

* **proximal + anode**: what the patient gets;
* **anode alone**: does the anode itself fire the fibre (virtual cathodes, anode break)? It tells *why*
  an AP reaches the muscle when it does.

Measured per setting, on the target fascicle:

* `frac_muscle`: fraction of the fascicle's fibres that deliver an AP to the muscle end (criterion C1:
  below `BLOCK_LIMIT`);
* `block_ratio`: of the fibres the proximal pulse fires, the fraction that no longer reach the muscle
  (the anodic block of Tovbis et al. 2026);
* `frac_anode_fires`: fraction of the fascicle fired by the anode alone;
* `frac_fired_at_fibre_end`: of those, the ones whose first AP started within 3 nodes of a fibre END.
  This is an **artefact check** (D-17): real fibres do not end there. In the first run of this notebook
  (grounded end planes, 60 mm fibres) the anode fired every F7 fibre at its distal end from 26 µA up,
  ~20 mm from the anode, which hid the block entirely. With insulated end planes and fibres running to
  1 mm from them this column should stay at 0; a warning is printed if it does not.

First the pulse shapes (normalised to their amplitude), with the arrival of the volley.
''')

code(r'''
# Shapes of the anodal pulses tried (one per rise time), with the proximal pulse and the volley arrival.
fig, ax = plt.subplots(figsize=(8, 2.8))
for r in RISE_TIMES_s:
    p = anode_pulses(1.0, r)
    tt = np.array([[on, off] for on, off, _ in p]).ravel()
    aa = np.array([[a, a] for _, _, a in p]).ravel()
    ax.plot((tt - T_PROX_ON_s) * 1e3, aa, label=f"ramps {r*1e3:.1f} ms")
ax.axvspan(0, PROX_WIDTH_s * 1e3, color="k", alpha=0.3, label="proximal pulse")
ax.axvspan(np.nanmin(arr) * 1e3, np.nanmax(arr) * 1e3, color="r", alpha=0.3, label="volley at the anode")
ax.set(xlabel="time from the proximal onset [ms]", ylabel="anodal current / amplitude",
       title=f"Anodal pulses on {ANODE_NAME}")
ax.legend(fontsize=7, ncol=2)
fig.savefig(OUT_DIR / "anode_pulse_shapes.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
### Static check: what the return contacts change in the activating function

Before running NEURON: the activating function (second difference of the extracellular potential
along each fibre, per µA of the distal channel; > 0 depolarises) of the **monopolar** anode and of the
configuration chosen in `ANODE_RETURN`, computed from the FEM profiles only (superposition).

**How to read this figure.** *Left*: median over the fibres of the target fascicle. Around the distal
site, the negative dip is the hyperpolarisation that should block the volley, flanked by the two
virtual cathodes. Near the dashed lines at ±25 mm (the ends of the insulating sleeve) the monopolar
anode has a broad positive bump: that is where it fired every fascicle in the earlier runs. With a
balanced pattern the bump should be gone. *Right*: the largest value among the fibres of all the OTHER
fascicles: the return pads, cathodes during the anodal pulse, show up as positive peaks at the distal
site; they must stay well below what fires those fibres. The table gives the numbers.
''')

code(r'''
# Static screen without NEURON: activating function at the nodes for the monopolar anode and for this
# run's pattern (with its return contacts), near the distal site and near the sleeve end, target vs
# other fascicles. Positive AF = depolarisation (where the anode itself could fire fibres).
def af_nodes(i, pattern, amp_A=1e-6):
    """Activating function at the nodes of fibre i for amp_A of the pattern [V/m^2]; NaN at the end nodes."""
    ve = np.interp(FIBRE.z_node_m, z_m, fibre_jobs.pattern_profile(int(i), pattern)) * amp_A
    a = np.full_like(ve, np.nan)
    a[1:-1] = (ve[:-2] - 2 * ve[1:-1] + ve[2:]) / FIBRE.internode_m ** 2
    return a


configs = {"monopolar anode " + "+".join(f"AS{c}" for c in ANODE_CONTACTS): ANODE_CONTACTS}
if RETURN_PATTERN:
    configs[f"this run: {ANODE_NAME}"] = ANODE_PATTERN
others_idx = np.flatnonzero(fascicle != TARGET)
near_edge = np.abs(np.abs(FIBRE.z_node_m) - SLEEVE_HALF_m) < 2.5e-3
near_site = np.abs(FIBRE.z_node_m - Z_DIST_m) < 3e-3
fig, axs = plt.subplots(1, 2, figsize=(13, 4))
AF_ROWS = []
for (lab, pat), col in zip(configs.items(), ["grey", "tab:red"]):
    A = np.array([af_nodes(i, pat) for i in range(n_fibres)])
    axs[0].plot(zn_mm, np.nanmedian(A[TGT], 0), ".-", color=col, label=lab)
    axs[1].plot(zn_mm, np.nanmax(A[others_idx], 0), ".-", color=col, label=lab)
    AF_ROWS.append({"configuration": lab,
                    "sleeve ends: max depolarisation, all fibres": np.nanmax(A[:, near_edge]),
                    f"F{TARGET}: hyperpolarisation under the anode (median)": np.nanmedian(np.nanmin(A[TGT][:, near_site], 1)),
                    f"F{TARGET}: virtual cathodes beside the anode (median)": np.nanmedian(np.nanmax(A[TGT][:, near_site], 1)),
                    "other fascicles: max depolarisation at the distal site": np.nanmax(A[others_idx][:, near_site])})
for ax in axs:
    for s in (-1, 1):
        ax.axvline(s * SLEEVE_HALF_m * 1e3, color="k", ls="--", lw=0.8)
    ax.axvline(Z_DIST_m * 1e3, color="grey", ls=":", lw=0.8)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xlabel("z [mm]  (dashed: sleeve ends; dotted: distal site)")
    ax.legend(fontsize=7)
axs[0].set(ylabel="AF [V/m^2 per uA]", title=f"F{TARGET} (target): median over its fibres")
axs[1].set(title="other fascicles: maximum over their fibres")
fig.tight_layout()
fig.savefig(OUT_DIR / "anode_af_check.png", dpi=130, bbox_inches="tight")
plt.show()
AF_TABLE = pd.DataFrame(AF_ROWS).set_index("configuration")
AF_TABLE.to_csv(OUT_DIR / "anode_af_check.csv")
AF_TABLE.round(1).T
''')

code(r'''
# Stage 1, target fascicle only: for every (rise time, anodal current) two runs per fibre, anode alone and
# proximal + anode; per setting the residual volley at the muscle (C1), the block ratio, the fibres fired
# by the anode, the brain-side volley (C0) and the Shannon charge check.
GRID = [(rise, amp) for rise in RISE_TIMES_s for amp in ANODE_AMPS_A]
jobs = []
for rise, amp in GRID:
    pul = anode_pulses(amp, rise)
    t_stop = t_stop_for(rect(T_PROX_ON_s, PROX_WIDTH_s, I_PROX), pul)
    jobs += [(int(i), [(ANODE_PATTERN, pul)], t_stop) for i in TGT]
    jobs += [(int(i), [prox_drive(I_PROX), (ANODE_PATTERN, pul)], t_stop) for i in TGT]
print(f"{len(GRID)} settings x {len(TGT)} fibres x 2 runs = {len(jobs)} simulations on {N_WORKERS} processes ...")
t0 = time.time()
out = runner.run_ap_times(jobs)
print(f"done in {time.time()-t0:.0f} s")

n_t, rows, STAGE1_AP = len(TGT), [], {}
fired_t = FIRED[TGT]
for s, (rise, amp) in enumerate(GRID):
    alone, comb = out[2 * s * n_t:(2 * s + 1) * n_t], out[(2 * s + 1) * n_t:(2 * s + 2) * n_t]
    STAGE1_AP[(rise, amp)] = (alone, comb)
    nm_c = volley.counts_at(comb, NODE_M)
    fires_alone = (volley.counts_at(alone, NODE_B) + volley.counts_at(alone, NODE_M)) > 0
    rows.append({"rise_ms": rise * 1e3, "amp_uA": amp * 1e6, "x_Iprox": amp / I_PROX,
                 "frac_muscle": float(np.mean(nm_c > 0)),
                 "frac_brain": float(np.mean(volley.counts_at(comb, NODE_B) > 0)),
                 "block_ratio": float(np.mean(nm_c[fired_t] == 0)),
                 "frac_anode_fires": float(np.mean(fires_alone)),
                 "frac_fired_at_fibre_end": float(np.mean(volley.started_at_end(alone, N_NODES))),
                 "C1": bool(np.mean(nm_c > 0) < BLOCK_LIMIT),
                 **anode_safety(anode_pulses(amp, rise))})
STAGE1 = pd.DataFrame(rows)
STAGE1.to_csv(OUT_DIR / "stage1_target_fascicle.csv", index=False)
print(f"settings meeting C1 on F{TARGET}: {STAGE1.C1.sum()} of {len(STAGE1)}")
if STAGE1.frac_fired_at_fibre_end.max() > 0:
    print("WARNING (D-17): in some settings the anode fires fibres AT THEIR CUT END (first AP within 3 nodes of a "
          "fibre end), a model artefact: check the end boundary condition and the fibre length before trusting them.")
STAGE1.round(3)
''')

md(r'''
### How to read this figure: the block curves of the target fascicle

x axis: anodal current (log). One line per ramp length.

* *Left*: fraction of the target fascicle that still delivers an AP to the muscle. The grey dashed line
  is the value without anode (the recruited fraction); the dotted line is the C1 limit. A block shows as
  the curve falling below the dotted line; if it rises again at high current, the anode has started to
  fire the fibres itself (compare with the right panel).
* *Middle*: block ratio (blocked / fired by the proximal pulse); Tovbis et al. 2026 reach ≈ 0.9 at best.
* *Right*: fraction of the fascicle fired by the anode alone (virtual cathodes / anode break).

Hollow markers: settings beyond the Shannon limit (k = 1.75) on the anode contact.
''')

code(r'''
# Block curves of stage 1 against the anodal current (filled dots: Shannon-safe settings).
fig, axs = plt.subplots(1, 3, figsize=(15, 4), sharex=True)
for rise, g in STAGE1.groupby("rise_ms", sort=False):
    for ax, col in zip(axs, ["frac_muscle", "block_ratio", "frac_anode_fires"]):
        (ln,) = ax.semilogx(g.amp_uA, g[col], "-", label=f"ramps {rise:.1f} ms")
        ok = g.shannon_ok.to_numpy()
        ax.plot(g.amp_uA[ok], g[col][ok], "o", color=ln.get_color())
        ax.plot(g.amp_uA[~ok], g[col][~ok], "o", mfc="none", color=ln.get_color())
axs[0].axhline(OUT_BASE.frac_muscle[TARGET], color="grey", ls="--", lw=0.8, label="no anode")
axs[0].axhline(BLOCK_LIMIT, color="k", ls=":", lw=0.8, label="C1 limit")
axs[0].set(ylabel=f"fraction of F{TARGET} reaching the muscle", title="Residual volley at the muscle (C1)")
axs[1].set(ylabel="blocked / fired", title="Block ratio (Tovbis et al. 2026)", ylim=(-0.02, 1.02))
axs[2].set(ylabel=f"fraction of F{TARGET}", title="Fired by the anode alone")
for ax in axs:
    ax.set_xlabel(f"anodal current on {ANODE_NAME} [uA]")
axs[0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT_DIR / "stage1_block_curves.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
**Settings to check on the whole nerve.** If some settings meet C1, the `N_CHECK` with the lowest
charge are taken (least current, least side effects, safest). If none does, the `N_CHECK` with the
smallest residual volley are taken anyway, so that section 8 still shows how far the anode gets.
''')

code(r'''
# Settings passed to the whole-nerve check: those meeting C1, lowest charge first (or, if none, the
# smallest residual volley).
passing = STAGE1[STAGE1.C1].sort_values("charge_nC")
if len(passing):
    CANDS = passing.head(N_CHECK)
    print(f"{len(passing)} settings meet C1; checking the {len(CANDS)} with the lowest charge")
else:
    CANDS = STAGE1.sort_values(["frac_muscle", "charge_nC"]).head(N_CHECK)
    print("NO setting meets C1 on the target fascicle; checking the ones with the smallest residual volley")
CANDS[["rise_ms", "amp_uA", "x_Iprox", "frac_muscle", "block_ratio", "frac_anode_fires", "charge_nC",
       "density_uC_cm2", "shannon_ok"]].round(3)
''')

md(r'''
## 7. The best settings on the whole nerve

The same two runs (proximal + anode, anode alone) on **every** fibre, for the candidate settings.
Per fascicle (`volley.fascicle_outcomes`, compared with the proximal pulse alone):

* `frac_muscle` – fraction of the fascicle delivering an AP to the muscle end (C1 for the target, C2
  for the others);
* `frac_extra_brain` – fraction sending *more* APs to the brain than with the proximal pulse alone
  (C3: new sensations caused by the anode);
* `block_ratio` – blocked / fired by the proximal pulse.

A setting **passes** when C1, C2 and C3 hold. Safety is reported next to it (it is not part of the pass
flag, so that the physics and the hardware limits stay separate in the discussion).
''')

code(r'''
# Stage 2, whole nerve: every fibre with the anode alone and with proximal + anode, for each selected
# setting; criteria C0-C3 per fascicle (volley.fascicle_outcomes) and PASS = all four.
CHECK, OUT_FASC, AP_CHECK = [], {}, {}
for _, c in CANDS.iterrows():
    rise, amp = c.rise_ms * 1e-3, c.amp_uA * 1e-6
    pul = anode_pulses(amp, rise)
    t_stop = t_stop_for(rect(T_PROX_ON_s, PROX_WIDTH_s, I_PROX), pul)
    t0 = time.time()
    res = runner.run_ap_times([(i, [(ANODE_PATTERN, pul)], t_stop) for i in range(n_fibres)]
                              + [(i, [prox_drive(I_PROX), (ANODE_PATTERN, pul)], t_stop) for i in range(n_fibres)])
    alone, comb = res[:n_fibres], res[n_fibres:]
    AP_CHECK[(rise, amp)] = (alone, comb)
    o = volley.fascicle_outcomes(comb, fascicle, NODE_B, NODE_M, baseline=AP_BASE)
    OUT_FASC[(rise, amp)] = o
    oth = o.drop(index=TARGET)
    row = {"rise_ms": rise * 1e3, "amp_uA": amp * 1e6, "x_Iprox": amp / I_PROX,
           "target_frac_muscle": o.frac_muscle[TARGET], "target_block_ratio": o.block_ratio[TARGET],
           "C0_sensation_kept": bool(o.frac_brain[TARGET] >= 0.9 * OUT_BASE.frac_brain[TARGET]),
           "C1_target_blocked": o.frac_muscle[TARGET] < BLOCK_LIMIT,
           "C2_others_max_frac_muscle": oth.frac_muscle.max(),
           "C2_ok": bool((oth.frac_muscle < SIDE_LIMIT).all()),
           "C3_others_max_frac_extra_brain": oth.frac_extra_brain.max(),
           "C3_ok": bool((oth.frac_extra_brain < SIDE_LIMIT).all()),
           "fibres_fired_by_anode_alone": int(np.sum((volley.counts_at(alone, NODE_B) + volley.counts_at(alone, NODE_M)) > 0)),
           "of_which_at_fibre_end": int(np.sum(volley.started_at_end(alone, N_NODES))),
           **anode_safety(pul)}
    row["PASS"] = row["C0_sensation_kept"] and row["C1_target_blocked"] and row["C2_ok"] and row["C3_ok"]
    CHECK.append(row)
    o.to_csv(OUT_DIR / f"outcomes_anode_r{rise*1e3:.1f}ms_{amp*1e6:.0f}uA.csv")
    print(f"ramps {rise*1e3:.1f} ms, {amp*1e6:.0f} uA: PASS={row['PASS']} (C0 {row['C0_sensation_kept']}, C1 {row['C1_target_blocked']}, "
          f"C2 {row['C2_ok']}, C3 {row['C3_ok']}); Shannon ok: {row['shannon_ok']}  [{time.time()-t0:.0f} s]")
CHECK = pd.DataFrame(CHECK)
CHECK.to_csv(OUT_DIR / "stage2_whole_nerve.csv", index=False)
CHECK.round(3)
''')

md(r'''
### How to read this figure: per-fascicle outcome of each setting

One panel per setting. For each fascicle, two pairs of bars: what reaches the **muscle end**
(proximal pulse alone = light, with the anode = dark) and what reaches the **brain end** (same
convention). The dotted line is the 10 % level. A good setting brings the target's dark muscle bar
below the line, leaves its brain bar unchanged (the sensation is kept), and leaves every other
fascicle below the line at both ends.
''')

code(r'''
# Per-fascicle outcome of each checked setting: fraction reaching the muscle / brain end, with and
# without the anode (dotted line = 10 % limit).
fig, axs = plt.subplots(len(OUT_FASC), 1, figsize=(11, 3.2 * len(OUT_FASC)), squeeze=False, sharex=True)
x = np.arange(len(FASCICLES))
for ax, ((rise, amp), o) in zip(axs[:, 0], OUT_FASC.items()):
    ax.bar(x - 0.3, OUT_BASE.frac_muscle.loc[FASCICLES], 0.2, color="#ef9a9a", label="muscle end, proximal alone")
    ax.bar(x - 0.1, o.frac_muscle.loc[FASCICLES], 0.2, color="#c62828", label="muscle end, with anode")
    ax.bar(x + 0.1, OUT_BASE.frac_brain.loc[FASCICLES], 0.2, color="#90caf9", label="brain end, proximal alone")
    ax.bar(x + 0.3, o.frac_brain.loc[FASCICLES], 0.2, color="#1565c0", label="brain end, with anode")
    ax.axhline(RECRUIT, color="k", ls=":", lw=0.8)
    ax.set(ylabel="fraction of the fascicle", ylim=(0, 1.02),
           title=f"anode {ANODE_NAME}: ramps {rise*1e3:.1f} ms, {amp*1e6:.0f} uA ({amp/I_PROX:.1f} x I_prox)")
axs[-1, 0].set_xticks(x, [f"F{k}" + ("*" if k == TARGET else "") for k in FASCICLES])
axs[-1, 0].set_xlabel(f"fascicle (* = target F{TARGET})")
axs[0, 0].legend(fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.35))
fig.tight_layout()
fig.savefig(OUT_DIR / "stage2_per_fascicle.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
### How to read this figure: who is blocked, on the cross-section

Each dot is a fibre line, coloured by its fate with the best setting: **blocked** (fired by the
proximal pulse, no AP at the muscle), **passes** (fired, its AP still reaches the muscle), **fired by
the anode** (not fired by the proximal pulse, but active with the anode), **quiet**. The gold square is
the anode contact projected on the section. Expect the blocked fibres close to the anode and the ones
that pass on the far side of the fascicle (the hyperpolarisation decays with distance).
''')

code(r'''
# Cross-section map of the best setting: which fibres are blocked, still pass, are fired by the anode, or
# stay quiet.
best_row = CHECK[CHECK.PASS].sort_values("charge_nC").head(1)
if len(best_row) == 0:
    best_row = CHECK.sort_values("target_frac_muscle").head(1)
    print("no setting passes all criteria: showing the one with the smallest residual volley in the target")
BEST = (float(best_row.rise_ms.iloc[0]) * 1e-3, float(best_row.amp_uA.iloc[0]) * 1e-6)
BEST = min(AP_CHECK, key=lambda s: abs(s[0] - BEST[0]) + abs(s[1] - BEST[1]))   # exact dict key
alone_b, comb_b = AP_CHECK[BEST]
nm_b = volley.counts_at(comb_b, NODE_M)
act_b = (volley.counts_at(comb_b, NODE_B) + nm_b) > 0
status = np.where(FIRED & (nm_b == 0), "blocked", np.where(FIRED, "passes", np.where(act_b, "fired by the anode", "quiet")))
pal = {"blocked": "#2e7d32", "passes": "#c62828", "fired by the anode": "#f9a825", "quiet": "#d0d0d0"}
fig, ax = plt.subplots(figsize=(6.5, 5.5))
ax.scatter(fibre_x_m * 1e3, fibre_y_m * 1e3, s=10, c=[pal[s] for s in status])
for k in FASCICLES:
    sel = fascicle == k
    ax.text(fibre_x_m[sel].mean() * 1e3, fibre_y_m[sel].mean() * 1e3, f"F{k}", fontsize=8, weight="bold",
            ha="center", va="center")
for c in ANODE_PATTERN:
    s = sites_all.loc[c]
    if s["type"] == 1:                                  # shank contact (nominal layout position)
        ax.plot(s.r_mm * np.cos(np.radians(s.theta_deg)), s.r_mm * np.sin(np.radians(s.theta_deg)), "s",
                ms=9, mfc="gold", mec="k")
ax.set_aspect("equal")
ax.set(xlabel="x [mm]", ylabel="y [mm]",
       title=f"ramps {BEST[0]*1e3:.1f} ms, {BEST[1]*1e6:.0f} uA on {ANODE_NAME}\n(contact at its nominal layout position)")
ax.legend(handles=[Patch(color=v, label=k) for k, v in pal.items()], fontsize=7, loc="lower right")
fig.savefig(OUT_DIR / "best_cross_section.png", dpi=130, bbox_inches="tight")
plt.show()
print({k: int(np.sum(status[TGT] == k)) for k in pal}, f"(fibres of F{TARGET})")
''')

md(r'''
## 8. The volley with the block

The same pictures as in section 5.5, now with the best setting of section 7 (the one with the lowest
charge among those that pass; if none passes, the one with the smallest residual volley).

**How to read this figure.** *Left and middle*: the volley of the target fascicle without and with the
anode (space × time, same colour scale). With the block, the band going up (towards the muscle) fades
or stops at the anode line; the band going down (towards the brain) is unchanged. *Right*: the **reach
profile** of the target fascicle, the fraction of its fibres in which an AP reached each node. Without
the anode it is flat at the recruited fraction; with the anode it drops on the far side of the distal
site. The value on the right end is what the muscle receives (C1: below the dotted line). The GIF
shows the two volleys travelling together, plus the activity the anode adds in the other fascicles.
''')

code(r'''
# The volley of the target fascicle with and without the block (volley maps), the other fascicles with
# the block, and how far along the nerve the APs travel (reach profile); animated GIFs are saved too.
pul_b = anode_pulses(BEST[1], BEST[0])
t0_frames = min(T_PROX_ON_s, start_s(pul_b)) - 0.1e-3
FRAMES_B = np.arange(t0_frames, max(FRAMES_BASE[-1], end_s(pul_b) + 0.3e-3), FRAME_STEP_s)
others_idx = np.flatnonzero(fascicle != TARGET)
vm_no = volley.volley_map(AP_BASE, TGT, N_NODES, FRAMES_B, VOLLEY_WINDOW_s, n_norm=len(TGT))
vm_bl = volley.volley_map(comb_b, TGT, N_NODES, FRAMES_B, VOLLEY_WINDOW_s, n_norm=len(TGT))
vm_ot = volley.volley_map(comb_b, others_idx, N_NODES, FRAMES_B, VOLLEY_WINDOW_s, n_norm=len(others_idx))
reach_no = volley.reach_profile(AP_BASE, TGT, N_NODES, n_norm=len(TGT))
reach_bl = volley.reach_profile(comb_b, TGT, N_NODES, n_norm=len(TGT))

fig, axs = plt.subplots(1, 3, figsize=(16, 4.4), gridspec_kw={"width_ratios": [1, 1, 0.8]})
vmax = max(0.05, vm_no.max())
for ax, vm, ttl in [(axs[0], vm_no, "no anode"), (axs[1], vm_bl, f"anode {BEST[1]*1e6:.0f} uA, ramps {BEST[0]*1e3:.1f} ms")]:
    im = ax.pcolormesh((FRAMES_B - T_PROX_ON_s) * 1e3, zn_mm, vm.T, cmap="viridis", vmin=0, vmax=vmax, shading="auto")
    for zz, _ in SITE_MARKS:
        ax.axhline(zz * 1e3, color="w", lw=0.8, ls="--")
    ax.set(xlabel="time from the proximal onset [ms]", title=f"Volley of F{TARGET}: {ttl}")
axs[0].set_ylabel("node position z [mm]  (muscle at the top)")
fig.colorbar(im, ax=axs[:2], label="fraction of the fascicle with an AP")
axs[2].plot(zn_mm, reach_no, color="grey", ls="--", label="no anode")
axs[2].plot(zn_mm, reach_bl, color=colours[f"F{TARGET}"], lw=2, label="with anode")
axs[2].axhline(BLOCK_LIMIT, color="k", ls=":", lw=0.8, label="C1 limit")
axs[2].axvline(FIBRE.z_node_m[NODE_M] * 1e3, color="r", lw=0.8)
axs[2].text(FIBRE.z_node_m[NODE_M] * 1e3, 1.0, "muscle\nend", color="r", fontsize=7, ha="right", va="top")
for zz, _ in SITE_MARKS:
    axs[2].axvline(zz * 1e3, color="grey", ls="--", lw=0.8)
axs[2].set(xlabel="z [mm]", ylabel=f"fraction of F{TARGET} reached by an AP", ylim=(0, 1.02), title="Reach profile")
axs[2].legend(fontsize=7)
fig.savefig(OUT_DIR / "volley_with_block.png", dpi=130, bbox_inches="tight")
plt.show()

gif = volley.animate_profiles(OUT_DIR / "volley_with_block.gif",
                              {f"F{TARGET}, no anode": vm_no, f"F{TARGET}, with anode": vm_bl,
                               "other fascicles, with anode": vm_ot},
                              FRAMES_B, FIBRE.z_node_m, ylabel="fraction of the fascicle(s) with an AP",
                              ylim=(0, max(0.15, vm_no.max()) * 1.1), markers=SITE_MARKS, t_ref_s=T_PROX_ON_s,
                              fps=GIF_FPS, colours={f"F{TARGET}, no anode": "grey",
                                                     f"F{TARGET}, with anode": colours[f"F{TARGET}"],
                                                     "other fascicles, with anode": "#f9a825"},
                              styles={f"F{TARGET}, no anode": "--"},
                              title=f"Volley of F{TARGET} with and without the anode ({BEST[1]*1e6:.0f} uA on {ANODE_NAME}, "
                                    f"ramps {BEST[0]*1e3:.1f} ms)")
display(Image(filename=str(gif)))
''')

md(r'''
## 9. Summary and limitations

`summary.json` collects the numbers to quote in the report: target, maximum selective current and the
recruited fraction, conduction velocity and arrival time, the stage-1 table, the whole-nerve check of
the candidates, and the best setting.

**Limitations to keep in mind when reading the results** (also in `docs/block.md` §7 and D-16):

* one fibre diameter (10 µm), straight fibres, a node of Ranvier aligned with the proximal contact (the
  best case for excitation); smaller fibres conduct more slowly (the anode window must be longer) and
  are blocked at different currents (Tovbis et al.: block depends strongly on diameter and distance);
* motor and sensory fibres are not distinguished: "reaching the muscle end" is counted for every fibre;
* the 10 % limit for the residual volley at the muscle is an assumption carried over from the sensory
  threshold criterion (section 1.3), to confirm with the tutor;
* superposition with insulating idle contacts (D-07); quasi-static tissue; no electrode–tissue
  interface; the Shannon criterion is fitted on macro-electrodes and only indicative for small
  contacts, and the CIC of the gold contacts is not known yet.
''')

code(r'''
# Summary of the run (parameters, thresholds, criteria, best setting) written to summary.json.
def _clean(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


summary = {
    "patient": PATIENT, "fem_source": fem_source, "fibre_diam_um": FIBRE_DIAM_um,
    "prox_pattern": PROX_PATTERN, "anode_pattern": ANODE_PATTERN, "anode_return": RETURN_PATTERN,
    "anode_tag": ANODE_TAG, "target_fascicle": TARGET,
    "th10_target_uA": TH10[TARGET] * 1e6, "selective": bool(SELECTIVE), "I_prox_uA": I_PROX * 1e6,
    "I_prox_nC": I_PROX * PROX_WIDTH_s * 1e9, "recruited_target": REC.recruited_at_I_prox[TARGET],
    "fired_fibres": int(FIRED.sum()), "fired_target": int(FIRED[TGT].sum()),
    "conduction_velocity_mps": CV_mps, "arrival_at_distal_ms": ARRIVAL_s * 1e3,
    "anode_plateau_ms": [PLATEAU_START_s * 1e3, PLATEAU_END_s * 1e3],
    "criteria": {"BLOCK_LIMIT": BLOCK_LIMIT, "SIDE_LIMIT": SIDE_LIMIT, "shannon_k": SHANNON_K},
    "stage1_settings_meeting_C1": int(STAGE1.C1.sum()), "stage1": STAGE1.to_dict(orient="records"),
    "stage2": CHECK.to_dict(orient="records"),
    "best_setting": {"rise_ms": BEST[0] * 1e3, "amp_uA": BEST[1] * 1e6,
                     "passes": bool(CHECK[(CHECK.rise_ms == BEST[0] * 1e3) & (CHECK.amp_uA == BEST[1] * 1e6)].PASS.any())},
}
(OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=1, default=_clean))
runner.close()
print("wrote", OUT_DIR / "summary.json", "and stopped the worker processes")
''')

md(r'''
## 10. Export to the 3D viewer (Volley playback)

Writes the runs of this notebook to `results/viewer/<patient>/volley_<name>.js` for the viewer's
*Volley · NEURON playback* section (README, *Volley playback*): the proximal pulse alone, and for every
setting checked in section 7 the proximal pulse + anode and the anode alone. Only the AP times are
written (a few hundred kB per run). Then reload `viewer/index.html` and press <kbd>V</kbd>.
''')

code(r'''
# Export the runs to the 3D viewer (viewer/index.html, Volley playback): results/viewer/<patient>/
# volley_<name>.js + volleys.js (src/block/volley.py export_viewer_volley).
VIEW_DIR = ROOT / "results" / "viewer" / PATIENT
prox_ch = {"name": f"proximal {PROX_NAME}", "kind": "cathode", "contacts": PROX_PATTERN,
           "pulses": rect(T_PROX_ON_s, PROX_WIDTH_s, I_PROX)}


def anode_ch(rise, amp):
    return {"name": f"anode {ANODE_NAME}", "kind": "anode", "contacts": ANODE_PATTERN, "pulses": anode_pulses(amp, rise)}


common = dict(fibre_x_m=fibre_x_m, fibre_y_m=fibre_y_m, fascicle=fascicle, z_nodes_m=FIBRE.z_node_m,
              node_brain=NODE_B, node_muscle=NODE_M, t_ref_s=T_PROX_ON_s, target=TARGET,
              z_sites_m=[Z_PROX_m, Z_DIST_m])
volley.export_viewer_volley(VIEW_DIR, f"{PROX_NAME}_no_anode", f"{PROX_NAME} {I_PROX*1e6:.0f} uA, no anode",
                            AP_BASE, channels=[prox_ch], note=f"proximal pulse alone, MRG {FIBRE_DIAM_um} um", **common)
for (rise, amp), (alone, comb) in AP_CHECK.items():
    tag = f"{ANODE_TAG}_{amp*1e6:.0f}uA_r{rise*1e3:.1f}ms"
    passed = bool(CHECK[(CHECK.rise_ms == rise * 1e3) & (CHECK.amp_uA == amp * 1e6)].PASS.any())
    volley.export_viewer_volley(VIEW_DIR, f"{PROX_NAME}_{tag}",
                                f"{PROX_NAME} + anode {ANODE_NAME} {amp*1e6:.0f} uA, ramps {rise*1e3:.1f} ms",
                                comb, channels=[prox_ch, anode_ch(rise, amp)],
                                note=f"criteria C1-C3 {'met' if passed else 'NOT met'} (D-16)", **common)
    volley.export_viewer_volley(VIEW_DIR, f"{tag}_alone", f"anode {ANODE_NAME} alone {amp*1e6:.0f} uA, ramps {rise*1e3:.1f} ms",
                                alone, channels=[anode_ch(rise, amp)], note="anode without the proximal pulse", **common)
print("exported:", sorted(p.name for p in VIEW_DIR.glob("volley_*.js")))
''')

# =====================================================================================================
nb = {"cells": [], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                "language_info": {"name": "python"}},
      "nbformat": 4, "nbformat_minor": 5}
for i, (kind, src) in enumerate(cells):
    lines = src.split("\n")
    source = [l + "\n" for l in lines[:-1]] + [lines[-1]]
    cell = {"cell_type": kind, "id": f"c{i:02d}", "metadata": {}, "source": source}
    if kind == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    nb["cells"].append(cell)
OUT.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print("wrote", OUT, len(cells), "cells")
