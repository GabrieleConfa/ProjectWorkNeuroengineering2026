"""Generates notebooks/block_model_check.ipynb."""
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


md(r'''
# Why the anode did not block: a textbook check of the fibre model, and the slow-ramp block

**Where we are (2026-11-27).** `block_anode_design.ipynb` designed the best anodal configurations for 1–6
contacts and checked them in NEURON on fascicle F7: **none blocked a single fibre**. Every design went
from "0 % blocked" straight to "the distal pattern fires all 22 F7 fibres" between two current steps
(`results/block/anode_design/check_target_fascicle.csv`). The single mirror anode had failed the same way
(`docs/block.md` §6).

**The question (user).** The tutor and the theory say that an anode hyperpolarises a region and stops the
APs. Is the NEURON model wrong, or the code?

**What the literature says.** Bhadra & Kilgore 2004 (*IEEE TNSRE* 12:313, doi:10.1109/TNSRE.2004.834205)
studied direct-current (DC) block in NEURON (squid, frog and mammalian fibre models) and in the frog
sciatic nerve. Their findings (abstract):

* DC block works with both polarities, but **depolarisation block needs less current than
  hyperpolarisation block**;
* the block can happen **under the electrode or at the adjacent virtual-electrode sites** (for an anode:
  the depolarised flanks, where the sodium channels become inactivated);
* the electrode fires the fibre when the current is switched on quickly (the "onset response"). The
  tutor's picture (the anode hyperpolarises and stops the APs) is right, but the current needed is higher
  than the current that fires the fibre, unless the current rises slowly.

Every pulse we have tried so far rose in ≤ 1 ms and lasted ≤ 1 ms: exactly the regime where the fibre is
fired before it is blocked. This notebook checks that explanation and tries the cure.

What it does

* **A. Textbook check of our fibre model (no FEM, ~2 min).** One MRG fibre, point electrodes in a
  homogeneous medium (the analytic potential V = I / (4 π σ r)). An AP is launched upstream; a DC
  electrode downstream is switched on either abruptly or with a slow ramp, as an anode and as a cathode.
  If our model is right it must reproduce Bhadra & Kilgore: block exists with both polarities, and the
  cathode blocks at a lower current. If it does, the model and the code can block, and our failure comes
  from the stimulation, not from NEURON.
* **B. Where and when did the round-4 designs fire the fibres? (~1 min).** For each design, at the first
  current that fired F7: the node where each AP started and its time: at switch-on, during the plateau
  (virtual cathode), or after switch-off (anode break).
* **C. Slow-ramp block in the real nerve (~20–30 min).** The anode current rises slowly (2 or 8 ms), holds
  while the volley passes, and falls slowly, so that the fibres adapt instead of firing. Checked on F7, then the
  best settings on the whole nerve with the criteria C0–C3 (D-16), plus the charge safety.

Outputs: `results/block/model_check/`. Runs on all cores but one.
''')

md(r'''
## 1. Setup
''')

code(r'''
# Imports, project paths, NEURON check, output folder.
import json, os, sys, time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd().resolve()
if not (ROOT / "src" / "block").exists():
    ROOT = ROOT.parent
assert (ROOT / "src" / "block" / "fibre_jobs.py").exists(), "open the notebook from the notebooks/ folder"
from neuron import h
sys.path.insert(0, str(ROOT / "src" / "block"))
import fibre_jobs, volley

PATIENT = "UA01"
MECH_DIR = ROOT / "results" / "mrg" / "model"
BLOCK_DIR = ROOT / "results" / "block"
H5_PATH = BLOCK_DIR / f"{PATIENT}_fibre_profiles.h5"
DESIGN_DIR = BLOCK_DIR / "anode_design"
OUT_DIR = BLOCK_DIR / "model_check"
OUT_DIR.mkdir(parents=True, exist_ok=True)
P = json.loads((BLOCK_DIR / f"{PATIENT}_block_params.json").read_text())
MRG = P["mrg"]
FIBRE_DIAM_um = MRG["fibreDiam_um"]
N_WORKERS = max(1, (os.cpu_count() or 2) - 1)
print("NEURON", h.nrnversion(), "| outputs:", OUT_DIR)
''')

md(r'''
## 2. Parameters

* `SIGMA_S_m`: conductivity of the homogeneous medium of part A: endoneurium, longitudinal, 0.571 S/m
  (`docs/parameters.md`). It only scales the currents; part A reports every current in units of the
  cathodic threshold at the same distance, so the conclusions do not depend on it.
* `D_LIST_m`: distances electrode–fibre of part A.
* `RAMPS_A_s`, `RAMPS_C_s`: ramp durations. 0 = abrupt switch-on (what we did so far).
* `HOLD_POST_s`: how long the plateau lasts after the proximal pulse (the volley passes the distal site
  within ~0.5 ms and the muscle detection node within ~0.8 ms).
* `C_CONFIGS`: the configurations of part C. Single contacts close to F7 on the muscle side (the
  selectivity ranking of the contacts for F7 is printed in section 5), plus two designs of
  `block_anode_design.ipynb`. Weights: sign = polarity (anode > 0), the largest |weight| carries `amp`.
''')

code(r'''
SIGMA_S_m = 0.571                     # endoneurium, longitudinal (docs/parameters.md); only scales currents
D_LIST_m = [0.25e-3, 0.5e-3, 1.0e-3]  # electrode-fibre distance, part A
Z_LAUNCH_m, D_LAUNCH_m = -20e-3, 0.5e-3   # part A: AP launched by a cathode upstream
Z_BLOCK_m = 5e-3                      # part A: the DC electrode, over the node closest to this z
RAMPS_A_s = [0.0, 5e-3]               # part A: abrupt vs slow switch-on
FACTORS_A = np.geomspace(0.3, 100, 18)    # part A: DC current / cathodic threshold (50 us) at the same distance
RAMP_STEP_s = 40e-6                   # the ramps are staircases of 40 us steps
T0_s = 0.2e-3                         # start of the distal current
HOLD_PRE_s, HOLD_POST_s = 0.3e-3, 1.5e-3
T_AFTER_s = P["pulse"]["tAfter_s"]
STIM_W_s = P["pulse"]["stimWidth_s"]  # proximal / launch pulse width (50 us)

C_CONFIGS = {                         # part C: name -> {contact: weight}
    "AS27": {27: 1.0}, "AS26": {26: 1.0}, "AS33": {33: 1.0}, "AS32": {32: 1.0}, "AS38": {38: 1.0},
}
C_DESIGNS = ["free_K2_N3", "free_K2_N4"]  # added from results/block/anode_design/candidates.json
POLARITY = +1                         # +1 anodal (the tutor's scheme); -1 tries the same contacts as cathodes
RAMPS_C_s = [2e-3, 8e-3]
AMPS_C_A = np.geomspace(20e-6, 2e-3, 10)
N_CHECK = 3                           # settings re-checked on the whole nerve
SHANNON_K = 1.75
CONTACT_AREA_m2 = 0.2e-3 * 0.2e-3
''')

md(r'''
## 3. Shared pieces

One MRG fibre (10 µm, 68 mm, a node at the proximal station z = −10.3 mm, as in the other notebooks) is
built in every worker process. `fibre_jobs.run_profiles` runs jobs whose potential along the fibre is given
explicitly (part A), `run_ap_times` jobs on the FEM profiles (parts B and C).

`ramp_hold` builds a current waveform: linear rise over `ramp`, plateau `hold`, linear fall over `ramp`,
written as a staircase of short constant steps (NEURON applies piecewise-constant currents).
''')

code(r'''
# Shared pieces: the MRG fibre of the nerve runs (node under AS 2), pulse helpers (ramp + hold written as
# a staircase), time of the first AP at a node, and the Shannon charge check on the most loaded contact.
with h5py.File(H5_PATH, "r") as f:
    CZ = dict(zip(f["contact_id"][:].ravel().astype(int).tolist(), f["contact_z_m"][:].ravel().tolist()))
    fascicle = f["fascicle"][:].ravel().astype(int)
n_fibres = len(fascicle)
runner = fibre_jobs.FibreRunner(H5_PATH, MRG, FIBRE_DIAM_um, CZ[2], MECH_DIR, N_WORKERS)
FIBRE = runner.fibre
ZN = FIBRE.z_node_m
NODE_B, NODE_M = runner.node_brain, runner.node_muscle
print(f"{FIBRE.n_nodes} nodes; detection nodes at z = {ZN[NODE_B]*1e3:.1f} (brain) and {ZN[NODE_M]*1e3:.1f} mm (muscle)")


def ramp_hold(t0, ramp, hold, amp):
    """Pulses [(t_on, t_off, amplitude)]: rise over `ramp`, plateau `hold`, fall over `ramp`."""
    p = []
    n = int(np.ceil(ramp / RAMP_STEP_s)) if ramp > 0 else 0
    if n:
        e = t0 + np.arange(n + 1) * ramp / n
        p += [(e[k], e[k + 1], amp * (k + 0.5) / n) for k in range(n)]
    t1 = t0 + ramp
    t2 = t1 + hold
    p.append((t1, t2, amp))
    if n:
        e = t2 + np.arange(n + 1) * ramp / n
        p += [(e[k], e[k + 1], amp * (1 - (k + 0.5) / n)) for k in range(n)]
    return p


def t_end(pulses):
    return max(off for _, off, _ in pulses)


def first_t(ap_fibre, node):
    t = ap_fibre[node]
    return float(t[0]) if len(t) else np.nan


def safety(weights, pulses):
    """Charge per phase on the contact with the largest |weight| and the Shannon check (k = 1.75)."""
    q = sum((off - on) * a for on, off, a in pulses) * max(abs(w) for w in weights.values())
    ok, Q_uC, D = volley.shannon_ok(q, CONTACT_AREA_m2, SHANNON_K)
    return {"charge_nC": Q_uC * 1e3, "density_uC_cm2": D, "shannon_ok": ok}
''')

md(r'''
## 4. Part A — textbook check of the fibre model

**Set-up.** The fibre lies along z in an infinite homogeneous medium. A point electrode at distance d from
the fibre axis, over a node at z_e, gives along the fibre V(z) = I / (4 π σ √((z − z_e)² + d²)) (for +1 A,
anodal; a cathode is the same with the opposite sign). Three electrodes:

1. a **launch cathode** 0.5 mm from the fibre at z ≈ −20 mm, 50 µs at 2× its threshold: it starts an AP
   travelling towards the muscle end (it plays the role of the proximal stimulation);
2. the **DC electrode** at z ≈ +5 mm, distance d: an anode or a cathode, switched on abruptly or with a 5 ms
   ramp, held while the AP passes, then switched off in the same way.

Per setting two runs: DC electrode alone ("does it fire the fibre?"), and launch + DC electrode ("does
the launched AP still reach the muscle detection node?"). The launched AP is recognised by its arrival
time (window around its arrival without the DC electrode).

**What our model must show (Bhadra & Kilgore 2004):** a current range where the launched AP is blocked
for both polarities, the cathodal block starting at a lower current than the anodal one; with an abrupt
switch-on the DC electrode fires the fibre (onset response), and a slow ramp should avoid it.
''')

code(r'''
ZS = np.linspace(-MRG["fibreLength_m"] / 2, MRG["fibreLength_m"] / 2, 6801)   # 10 um samples


def point_V(z_e, d):
    """Potential [V] along the fibre for +1 A from a point electrode (homogeneous medium)."""
    return 1.0 / (4 * np.pi * SIGMA_S_m * np.sqrt((ZS - z_e) ** 2 + d ** 2))


zL = ZN[np.argmin(np.abs(ZN - Z_LAUNCH_m))]
zB = ZN[np.argmin(np.abs(ZN - Z_BLOCK_m))]
VL = -point_V(zL, D_LAUNCH_m)                 # launch electrode, cathodic


def sweep_threshold(profile, t_on, width, amps):
    """Smallest amplitude of `amps` whose single pulse fires the fibre (NaN if none)."""
    jobs = [([(ZS, profile, [(t_on, t_on + width, a)])], t_on + width + T_AFTER_s) for a in amps]
    out = runner.run_profiles(jobs)
    fired = [sum(len(t) for t in ap) > 0 for ap in out]
    return float(amps[fired.index(True)]) if any(fired) else np.nan


AMPS_TH = np.geomspace(1e-6, 20e-3, 80)       # 13 % steps
t0 = time.time()
I_LAUNCH = 2 * sweep_threshold(VL, 0.5e-3, STIM_W_s, AMPS_TH)
I_C50 = {d: sweep_threshold(-point_V(zB, d), 0.5e-3, STIM_W_s, AMPS_TH) for d in D_LIST_m}
print(f"launch current {I_LAUNCH*1e6:.0f} uA; cathodic thresholds (50 us) at the DC site: "
      + ", ".join(f"d = {d*1e3:.2f} mm: {I*1e6:.0f} uA" for d, I in I_C50.items()) + f"  [{time.time()-t0:.0f} s]")

# arrival of the launched AP without the DC electrode, and the conduction velocity of our fibre
T_L_REF = 0.5e-3
ref = runner.run_profiles([([(ZS, VL, [(T_L_REF, T_L_REF + STIM_W_s, I_LAUNCH)])], T_L_REF + 4e-3)])[0]
DELAY_M = first_t(ref, NODE_M) - T_L_REF
i1, i2 = np.argmin(np.abs(ZN - 0)), np.argmin(np.abs(ZN - 20e-3))
CV = (ZN[i2] - ZN[i1]) / (first_t(ref, i2) - first_t(ref, i1))
print(f"launched AP reaches the muscle node {DELAY_M*1e3:.2f} ms after the launch; conduction velocity {CV:.1f} m/s")
''')

code(r'''
rowsA, jobs, keys = [], [], []
for pol in (+1, -1):
    for d in D_LIST_m:
        VB = pol * point_V(zB, d)
        for ramp in RAMPS_A_s:
            for fac in FACTORS_A:
                amp = fac * I_C50[d]
                pul = ramp_hold(T0_s, ramp, HOLD_PRE_s + HOLD_POST_s + 0.5e-3, amp)
                t_launch = T0_s + ramp + HOLD_PRE_s
                ts = t_end(pul) + T_AFTER_s
                jobs.append(([(ZS, VB, pul)], ts))                                            # DC alone
                jobs.append(([(ZS, VB, pul), (ZS, VL, [(t_launch, t_launch + STIM_W_s, I_LAUNCH)])], ts))  # + launch
                keys.append((pol, d, ramp, fac, t_launch))
print(f"{len(jobs)} runs ...")
t0 = time.time()
outA = runner.run_profiles(jobs)
print(f"done in {time.time()-t0:.0f} s")
for s, (pol, d, ramp, fac, t_launch) in enumerate(keys):
    alone, comb = outA[2 * s], outA[2 * s + 1]
    t_arr = t_launch + DELAY_M
    tm = comb[NODE_M]
    passed = bool(np.any((tm > t_arr - 0.3e-3) & (tm < t_arr + 1.0e-3)))
    fired = sum(len(t) for t in alone) > 0
    t_fire = min((t[0] for t in alone if len(t)), default=np.nan)
    rowsA.append({"polarity": "anode" if pol > 0 else "cathode", "d_mm": d * 1e3, "ramp_ms": ramp * 1e3,
                  "factor": fac, "amp_uA": fac * I_C50[d] * 1e6, "launched_AP_blocked": not passed,
                  "fired_by_DC": fired, "first_AP_by_DC_ms": t_fire * 1e3,
                  "clean_block": (not passed) and (not fired)})
A_TAB = pd.DataFrame(rowsA)
A_TAB.to_csv(OUT_DIR / "partA_point_electrode.csv", index=False)


def first_factor(g, col):
    v = g[g[col]].factor
    return float(v.min()) if len(v) else np.nan


A_SUM = (A_TAB.groupby(["polarity", "d_mm", "ramp_ms"])
         .apply(lambda g: pd.Series({"block_from_x": first_factor(g, "launched_AP_blocked"),
                                     "fires_from_x": first_factor(g, "fired_by_DC"),
                                     "clean_block_from_x": first_factor(g, "clean_block"),
                                     "clean_block_settings": int(g.clean_block.sum())}))
         .reset_index())
A_SUM.to_csv(OUT_DIR / "partA_summary.csv", index=False)
A_SUM.round(2)
''')

md(r'''
**How to read the table.** Currents are in units of the cathodic threshold (50 µs pulse) at the same
distance. `block_from_x`: lowest current at which the launched AP no longer reaches the muscle end.
`fires_from_x`: lowest current at which the DC electrode alone fires the fibre. `clean_block_from_x`:
lowest current that blocks without firing (NaN: never). Check against Bhadra & Kilgore 2004: block for both
polarities, cathode before anode, and firing with the abrupt switch-on.

**How to read the figure.** One panel per polarity and ramp; x = current / cathodic threshold; solid =
launched AP blocked (1) or not (0), dashed = the DC electrode alone fires the fibre; one colour per
distance. A clean block is where the solid line is up and the dashed line is down.
''')

code(r'''
# Part A figure: for each polarity and ramp, whether the launched AP was blocked (solid) and whether the
# DC pulse itself fired the fibre (dashed), against the current factor, one line per electrode distance.
fig, axs = plt.subplots(2, len(RAMPS_A_s), figsize=(5.5 * len(RAMPS_A_s), 6.5), sharex=True, sharey=True, squeeze=False)
for i, pol in enumerate(["anode", "cathode"]):
    for j, ramp in enumerate(RAMPS_A_s):
        ax = axs[i, j]
        for d in D_LIST_m:
            g = A_TAB[(A_TAB.polarity == pol) & (A_TAB.ramp_ms == ramp * 1e3) & (A_TAB.d_mm == d * 1e3)]
            (ln,) = ax.semilogx(g.factor, g.launched_AP_blocked.astype(float) + 0.02 * D_LIST_m.index(d), "-o", ms=3,
                                label=f"d = {d*1e3:.2f} mm")
            ax.semilogx(g.factor, g.fired_by_DC.astype(float) - 0.02 * (1 + D_LIST_m.index(d)), "--", color=ln.get_color())
        ax.set_title(f"{pol}, ramp {ramp*1e3:.0f} ms", fontsize=9)
        ax.set_yticks([0, 1]); ax.set_yticklabels(["no", "yes"])
axs[1, 0].set_xlabel("DC current / cathodic threshold (50 us)")
axs[0, 0].legend(fontsize=7)
fig.suptitle("Part A: solid = launched AP blocked, dashed = DC electrode fires the fibre", fontsize=10)
fig.tight_layout()
fig.savefig(OUT_DIR / "partA_point_electrode.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
## 5. Part B — where and when did the round-4 designs fire the fibres?

For each design of `block_anode_design.ipynb`, at the first current of its sweep that fired F7: the design
alone (rectangular pulse, 0.8 ms, as in that notebook), on every F7 fibre. For each fibre: the node where
its first AP started and the time of that AP relative to the pulse:

* **onset** (within 0.2 ms of switch-on): the fast switch-on fired it;
* **plateau**: a depolarised region (virtual cathode, sleeve end) reached threshold during the pulse;
* **offset** (after switch-off): anode break: the hyperpolarised nodes rebound and fire.

Also printed: the ranking of the contacts by how selectively their cathodic activating function reaches F7
(used to pick the contacts of part C).
''')

code(r'''
CANDS = json.loads((DESIGN_DIR / "candidates.json").read_text())
SUM_D = json.loads((DESIGN_DIR / "summary.json").read_text())
TARGET, I_PROX = int(SUM_D["target"]), SUM_D["I_prox_uA"] * 1e-6
S1 = pd.read_csv(DESIGN_DIR / "check_target_fascicle.csv")
TGT_ALL = np.flatnonzero(fascicle == TARGET)
T_ON_D, W_D = 1.5e-3 + 0.1e-3, 0.8e-3            # pulse of block_anode_design.ipynb
rowsB, jobs, keys = [], [], []
for name, c in CANDS.items():
    g = S1[(S1.name == name) & (S1.frac_fired_by_pattern > 0)]
    if g.empty:
        continue
    amp = g.amp_uA.min() * 1e-6
    w = {int(k): v for k, v in c["weights"].items()}
    jobs += [(int(i), [(w, [(T_ON_D, T_ON_D + W_D, amp)])], T_ON_D + W_D + T_AFTER_s) for i in TGT_ALL]
    keys.append((name, amp, c["block_nodes"]))
t0 = time.time()
outB = runner.run_ap_times(jobs)
print(f"{len(jobs)} runs in {time.time()-t0:.0f} s")
nT = len(TGT_ALL)
for s, (name, amp, bn) in enumerate(keys):
    ap = outB[s * nT:(s + 1) * nT]
    org = volley.origin_nodes(ap)
    t1 = np.array([min((t[0] for t in a if len(t)), default=np.nan) for a in ap])
    ok = org >= 0
    phase = np.where(t1 < T_ON_D + 0.2e-3, "onset", np.where(t1 < T_ON_D + W_D, "plateau", "offset"))[ok]
    zo = ZN[org[ok]] * 1e3
    rowsB.append({"design": name, "amp_uA": amp * 1e6, "fired": int(ok.sum()), "of": nT,
                  "onset": int(np.sum(phase == "onset")), "plateau": int(np.sum(phase == "plateau")),
                  "offset": int(np.sum(phase == "offset")),
                  "origin_z_mm_median": float(np.median(zo)) if ok.any() else np.nan,
                  "origin_z_mm_range": f"{zo.min():.1f} .. {zo.max():.1f}" if ok.any() else "",
                  "block_nodes_z_mm": ", ".join(f"{ZN[k]*1e3:.1f}" for k in bn)})
B_TAB = pd.DataFrame(rowsB)
B_TAB.to_csv(OUT_DIR / "partB_where_when.csv", index=False)
B_TAB
''')

code(r'''
# selectivity of each contact for the target fascicle: 10%-value of the cathodic activating function (F7 / best other)
import anode_design
U = anode_design.unit_af(H5_PATH, ZN)
A_U, CID = U["A"], U["contact_id"]


def af10(a, k):
    v = np.sort(a[fascicle == k])[::-1]
    return v[max(0, int(np.ceil(0.1 * len(v))) - 1)]


rank = []
for r, c in enumerate(CID):
    a = (-A_U[r]).max(1)
    v = {k: af10(a, k) for k in np.unique(fascicle)}
    rank.append({"contact": int(c), "z_mm": CZ[int(c)] * 1e3, "sel_target": v[TARGET] / max(x for k, x in v.items() if k != TARGET)})
RANK = pd.DataFrame(rank).sort_values("sel_target", ascending=False)
RANK.head(12).round(2)
''')

md(r'''
## 6. Part C — slow-ramp anodal block in the nerve

The proximal pulse is the maximum selective pulse on AS 2 of the other notebooks (fires 13 of the 22 F7
fibres). The distal configuration ramps up over `ramp`, holds from 0.3 ms before the proximal pulse until
`HOLD_POST_s` after it, and ramps down. Two runs per F7 fibre: distal alone, distal + proximal.

Per setting (strict D-16 counting: ANY AP at the muscle node during the whole run counts):

* `frac_muscle` (C1: < 10 % of F7), `frac_brain` (C0: ≥ 90 % of the reference), `frac_fired_by_pattern`
  (the distal configuration alone fires the fibre), and `block_ratio_volley` (fraction of the fibres fired
  by the proximal pulse whose proximal AP does not reach the muscle on time: the block itself, even if the
  ramps fire the fibre at another moment);
* charge and Shannon check on the contact with the largest current. A slow ramp means a long phase, and
  the charge grows with it: this is the price of this approach and is reported honestly.
''')

code(r'''
for nm in C_DESIGNS:
    if nm in CANDS:
        C_CONFIGS[nm] = {int(k): v for k, v in CANDS[nm]["weights"].items()}
C_CONFIGS = {k: {c: POLARITY * w for c, w in v.items()} for k, v in C_CONFIGS.items()}
print("configurations:", C_CONFIGS)

T_BASE = 1.5e-3
prox_base = ({2: -1.0}, [(T_BASE, T_BASE + STIM_W_s, I_PROX)])
t0 = time.time()
AP_BASE = runner.run_ap_times([(i, [prox_base], T_BASE + STIM_W_s + T_AFTER_s) for i in range(n_fibres)])
OUT_BASE = volley.fascicle_outcomes(AP_BASE, fascicle, NODE_B, NODE_M)
FIRED = (volley.counts_at(AP_BASE, NODE_B) + volley.counts_at(AP_BASE, NODE_M)) > 0
DELAY = np.array([first_t(a, NODE_M) - T_BASE for a in AP_BASE])     # proximal pulse -> muscle node
print(f"reference: {FIRED[TGT_ALL].sum()} of {len(TGT_ALL)} F{TARGET} fibres fire  [{time.time()-t0:.0f} s]")


def drive_c(w, ramp, amp):
    pul = ramp_hold(T0_s, ramp, HOLD_PRE_s + HOLD_POST_s, amp)
    t_prox = T0_s + ramp + HOLD_PRE_s
    return pul, t_prox, t_end(pul) + T_AFTER_s


jobs, keys = [], []
for name, w in C_CONFIGS.items():
    for ramp in RAMPS_C_s:
        for amp in AMPS_C_A:
            pul, t_prox, ts = drive_c(w, ramp, amp)
            prox = ({2: -1.0}, [(t_prox, t_prox + STIM_W_s, I_PROX)])
            jobs += [(int(i), [(w, pul)], ts) for i in TGT_ALL] + [(int(i), [prox, (w, pul)], ts) for i in TGT_ALL]
            keys.append((name, ramp, amp, t_prox))
print(f"{len(jobs)} simulations on {N_WORKERS} processes (longer runs than before: ~{2*max(RAMPS_C_s)*1e3+5:.0f} ms each) ...")
t0 = time.time()
outC = runner.run_ap_times(jobs)
print(f"done in {time.time()-t0:.0f} s")
''')

code(r'''
# Part C table: per configuration, ramp and current, the residual volley at the muscle, the brain-side
# volley, the fibres fired by the pattern and the charge; criteria C0 and C1 on the target fascicle.
fired_t = FIRED[TGT_ALL]
rowsC = []
for s, (name, ramp, amp, t_prox) in enumerate(keys):
    alone, comb = outC[2 * s * nT:(2 * s + 1) * nT], outC[(2 * s + 1) * nT:(2 * s + 2) * nT]
    nm, nb = volley.counts_at(comb, NODE_M), volley.counts_at(comb, NODE_B)
    on_time = np.array([np.any(np.abs(comb[k][NODE_M] - (t_prox + DELAY[i])) < 0.5e-3) if np.isfinite(DELAY[i]) else False
                        for k, i in enumerate(TGT_ALL)])
    w = C_CONFIGS[name]
    rowsC.append({"name": name, "ramp_ms": ramp * 1e3, "amp_uA": amp * 1e6,
                  "frac_muscle": float(np.mean(nm > 0)), "frac_brain": float(np.mean(nb > 0)),
                  "frac_fired_by_pattern": float(np.mean((volley.counts_at(alone, NODE_B) + volley.counts_at(alone, NODE_M)) > 0)),
                  "block_ratio_volley": float(np.mean(~on_time[fired_t])),
                  **safety(w, drive_c(w, ramp, amp)[0])})
C_TAB = pd.DataFrame(rowsC)
C_TAB["C0"] = C_TAB.frac_brain >= 0.9 * OUT_BASE.frac_brain[TARGET]
C_TAB["C1"] = C_TAB.frac_muscle < P["recruitFraction"]
C_TAB.to_csv(OUT_DIR / "partC_target_fascicle.csv", index=False)
print("settings meeting C0 and C1 on F7:", int((C_TAB.C0 & C_TAB.C1).sum()))
C_TAB.round(3)
''')

md(r'''
**How to read the figure.** One panel per configuration, one colour per ramp. *Solid*: fraction of F7
reaching the muscle (any AP, C1 below the dotted line). *Dotted with markers*: block of the volley itself
(fraction of the fired fibres whose proximal AP does not arrive). *Dashed*: fraction fired by the distal
configuration alone. A working setting: solid below 10 %, dashed at 0, dotted near 1.
''')

code(r'''
# Part C figure: residual volley (solid), block ratio (dotted) and fibres fired by the pattern alone
# (dashed) against the current, one panel per configuration.
fig, axs = plt.subplots(1, len(C_CONFIGS), figsize=(3.8 * len(C_CONFIGS), 3.6), sharey=True, squeeze=False)
for ax, name in zip(axs[0], C_CONFIGS):
    g = C_TAB[C_TAB.name == name]
    for ramp, gr in g.groupby("ramp_ms"):
        (ln,) = ax.semilogx(gr.amp_uA, gr.frac_muscle, "-", label=f"ramp {ramp:.0f} ms")
        ax.semilogx(gr.amp_uA, gr.block_ratio_volley, ":o", ms=3, color=ln.get_color())
        ax.semilogx(gr.amp_uA, gr.frac_fired_by_pattern, "--", color=ln.get_color())
    ax.axhline(P["recruitFraction"], color="k", ls=":", lw=0.8)
    ax.set_title(name, fontsize=9)
    ax.set_xlabel("current of the largest contact [uA]")
axs[0, 0].set_ylabel(f"fraction of F{TARGET}")
axs[0, 0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT_DIR / "partC_target_curves.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
### 6.1 Whole-nerve check of the best settings

The `N_CHECK` settings meeting C0 and C1 on F7, lowest charge first; every fibre of the nerve, distal alone
and with the proximal pulse; pass = C0–C3 (D-16).
''')

code(r'''
# Whole-nerve check (C0-C3) of the best Part C settings.
SEL = C_TAB[C_TAB.C0 & C_TAB.C1].sort_values("charge_nC").head(N_CHECK)
if SEL.empty:
    print("no setting meets C0 and C1 on F7: nothing to check on the whole nerve")
CHECK = []
for _, s in SEL.iterrows():
    w = C_CONFIGS[s["name"]]
    pul, t_prox, ts = drive_c(w, s.ramp_ms * 1e-3, s.amp_uA * 1e-6)
    prox = ({2: -1.0}, [(t_prox, t_prox + STIM_W_s, I_PROX)])
    t0 = time.time()
    res = runner.run_ap_times([(i, [(w, pul)], ts) for i in range(n_fibres)] + [(i, [prox, (w, pul)], ts) for i in range(n_fibres)])
    comb = res[n_fibres:]
    o = volley.fascicle_outcomes(comb, fascicle, NODE_B, NODE_M, baseline=AP_BASE)
    oth = o.drop(index=TARGET)
    row = {"name": s["name"], "ramp_ms": s.ramp_ms, "amp_uA": s.amp_uA,
           "C0": bool(o.frac_brain[TARGET] >= 0.9 * OUT_BASE.frac_brain[TARGET]),
           "C1": bool(o.frac_muscle[TARGET] < P["recruitFraction"]),
           "C2": bool((oth.frac_muscle < P["recruitFraction"]).all()),
           "C3": bool((oth.frac_extra_brain < P["recruitFraction"]).all()),
           "worst_other_muscle": float(oth.frac_muscle.max()), "worst_other_extra_brain": float(oth.frac_extra_brain.max()),
           "charge_nC": s.charge_nC, "density_uC_cm2": s.density_uC_cm2, "shannon_ok": s.shannon_ok}
    row["PASS"] = row["C0"] and row["C1"] and row["C2"] and row["C3"]
    CHECK.append(row)
    o.to_csv(OUT_DIR / f"outcomes_{s['name']}_ramp{s.ramp_ms:.0f}ms_{s.amp_uA:.0f}uA.csv")
    print(f"{s['name']} ramp {s.ramp_ms:.0f} ms {s.amp_uA:.0f} uA: PASS={row['PASS']}  [{time.time()-t0:.0f} s]")
CHECK = pd.DataFrame(CHECK)
CHECK.to_csv(OUT_DIR / "partC_whole_nerve.csv", index=False)
CHECK
''')

md(r'''
## 7. Summary

Written to `results/block/model_check/summary.json`. The next step depends on
which of these holds.

* Part A shows a block window → the model and the code can block; the round-4 failure was the fast
  switch-on.
* Part B says onset/plateau/offset → which mechanism fired the fibres (and therefore what to fix).
* Part C passes → a working block in the nerve, with its charge cost.
''')

code(r'''
# Summary of parts A-C written to summary.json and printed.
summary = {"partA": A_SUM.to_dict(orient="records"), "cv_m_s": CV, "partB": B_TAB.to_dict(orient="records"),
           "partC_best_target": C_TAB[C_TAB.C0 & C_TAB.C1].sort_values("charge_nC").head(5).to_dict(orient="records"),
           "partC_whole_nerve": CHECK.to_dict(orient="records") if len(CHECK) else []}
(OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
ok_A = A_SUM.groupby("polarity").clean_block_settings.sum()
print("A: settings with a clean block:", ok_A.to_dict())
print("B: fibres fired at onset / plateau / offset:", B_TAB[["onset", "plateau", "offset"]].sum().to_dict())
if len(CHECK) and CHECK.PASS.any():
    b = CHECK[CHECK.PASS].iloc[0]
    print(f"C: PASS with {b['name']}, ramp {b.ramp_ms:.0f} ms, {b.amp_uA:.0f} uA, {b.charge_nC:.0f} nC/phase, "
          f"{b.density_uC_cm2:.0f} uC/cm^2 (Shannon ok: {b.shannon_ok})")
else:
    print("C: no setting passes C0-C3 yet")
runner.close()
''')

nb = {"cells": [], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
for i, (kind, src) in enumerate(cells):
    lines = src.split("\n")
    cell = {"cell_type": kind, "id": f"m{i:02d}", "metadata": {}, "source": [l + "\n" for l in lines[:-1]] + [lines[-1]]}
    if kind == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    nb["cells"].append(cell)
OUT.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print("wrote", OUT, len(cells), "cells")
