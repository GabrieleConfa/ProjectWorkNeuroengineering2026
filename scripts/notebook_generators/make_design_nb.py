"""Generates notebooks/block_anode_design.ipynb."""
import json
import sys
from pathlib import Path

OUT = Path(sys.argv[1])
cells = []


def md(s):
    cells.append(("markdown", s.strip("\n")))


def code(s):
    cells.append(("code", s.strip("\n")))


md(r'''
# Designing the distal block: the fewest contacts that stop the volley

**Question (user, 2026-11-27).** The mirror anode (the same contact as the proximal cathode, reversed)
never blocked the volley of the target fascicle (`docs/block.md` §6, third round). Abandon the mirror
rule and use any contact of the electrode: **which configuration blocks the volley with the fewest active
contacts and the least current?**

What this notebook does

1. **Design (static, seconds).** For every number of contacts N = 1, 2, 3, ... it finds, by optimisation
   on the FEM solution, the configuration (which contacts, signed currents) that best does what a block
   needs (section 1), and draws the trade-off "number of contacts vs block quality". No NEURON here.
2. **Check (NEURON, MRG fibres, ~10–15 min).** Each design is driven with increasing current together with
   the proximal pulse, on every fibre of the target fascicle; then the best settings on the whole nerve,
   with the fascicle-level criteria C0–C3 (D-16, C0 added 2026-11-27) and the charge limits.
3. **Result.** The smallest configuration that passes, written to
   `results/block/anode_design/candidates.json`; `notebooks/block_volley.ipynb` can then replay it with all
   the pictures (`DIST_PATTERN = "design:<name>"`).

Before running: the fibre profiles of the current solution (`export_block_profiles`) and the threshold
cache of `block_volley.ipynb` (`results/block/volley/UA01_AS2_D10_thresholds.csv`; if missing it is
recomputed, ~4 min). Outputs: `results/block/anode_design/`.
''')

md(r'''
## 1. What a block needs, written as an optimisation problem

From the earlier rounds (`docs/block.md` §6):

* **Several consecutive nodes.** A myelinated fibre conducts with a safety margin: an AP jumps over a single
  node that cannot fire. A 0.2 mm contact hyperpolarises about one node (internode 1.15 mm for 10 µm
  fibres), so the block must hyperpolarise **K consecutive nodes** of every fibre that the proximal pulse
  fires in the target fascicle.
* **No fibre may be fired by the distal site.** Wherever the configuration depolarises a fibre (virtual
  cathodes beside an anode, cathodic return contacts, the ends of the insulating sleeve), it must stay
  below what fires that fibre, in **every fascicle**.
* **The sensation must pass.** On the brain side of the stimulation site the target fibres must not be
  hyperpolarised (here: at most 30 % of the block depth).

The model is linear: the activating function (AF, second difference of the extracellular potential along
the fibre, > 0 depolarises) of a configuration with currents I is AF = A · I, where A holds the AF of every
contact for 1 A (one FEM solution per contact). So the problem is linear in I:

*maximise h* such that  AF ≤ 1 at every node of every fibre,  AF ≤ −h at the K block nodes of every
target fibre,  AF ≥ −0.3 h on their brain side,  and at most N contacts carry current.

h is then the ratio **block depth / worst depolarisation**: it does not depend on the (unknown) absolute AF
thresholds, only on the shape of the field. The "at most N contacts" condition needs a yes/no variable per
contact: a mixed-integer linear program, solved with `scipy.optimize.milp` (HiGHS). The optimisation also
chooses where along the nerve (between the two sites) to put the K block nodes.

Limits of the design step: the AF is a proxy. It ignores the pulse shape and the width of a depolarised
region (a weak but broad bump can fire a fibre under a long pulse: the sleeve-end effect), and the actual
ratio needed for a block is not known. That is why every design is checked with the MRG fibres.
''')

md(r'''
## 2. Setup
''')

code(r'''
import json, math, os, sys, time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd().resolve()
if not (ROOT / "src" / "block").exists():
    ROOT = ROOT.parent
assert (ROOT / "src" / "block" / "anode_design.py").exists(), "open the notebook from the notebooks/ folder"
try:
    from neuron import h
except ImportError as err:
    raise ImportError("NEURON is not importable from this kernel (docs/block.md §4)") from err
sys.path.insert(0, str(ROOT / "src" / "block"))
import fibre_jobs, volley, anode_design

MECH_DIR = ROOT / "results" / "mrg" / "model"
BLOCK_DIR = ROOT / "results" / "block"
OUT_DIR = BLOCK_DIR / "anode_design"
OUT_DIR.mkdir(parents=True, exist_ok=True)
print("NEURON", h.nrnversion(), "| outputs:", OUT_DIR)
''')

md(r'''
## 3. Parameters

* `PROX_PATTERN`, `TARGET_FASCICLE`, `SEL_MARGIN`: as in `block_volley.ipynb` (maximum selective proximal
  pulse on AS 2; target chosen automatically).
* `K_LIST`: numbers of consecutive block nodes tried; `N_LIST`: maximum numbers of active contacts.
* `BLOCK_ZONE_m`: where the block nodes may be: distal of the proximal site and inside the sleeve.
* `EXCLUDE`: contacts the design may not use (the proximal stimulating contact).
* `DESIGN_VARIANTS`: "free" lets the net current go to the far ground or not, as the optimum prefers;
  "balanced" forces zero net current (optional: never better than "free").
* Pulse: the anodal plateau covers the arrival of the volley (0.1–0.9 ms after the proximal onset);
  `RISE_TIMES_s` ramps; `AMPS_A` = current on the contact with the largest weight.
''')

code(r'''
PATIENT = "UA01"
P = json.loads((BLOCK_DIR / f"{PATIENT}_block_params.json").read_text())
MRG = P["mrg"]
PROX_PATTERN = {2: -1.0}
TARGET_FASCICLE = None
SEL_MARGIN = 0.95
FIBRE_DIAM_um = MRG["fibreDiam_um"]
RECRUIT = P["recruitFraction"]

K_LIST = [1, 2, 3]                     # consecutive nodes to hyperpolarise
N_LIST = [1, 2, 3, 4, 6]               # maximum active contacts
BLOCK_ZONE_m = (-7e-3, 22e-3)          # block nodes between 3 mm distal of the proximal site and 3 mm
                                       # inside the sleeve end (z = 25 mm)
PROTECT_FRAC = 0.3                     # brain side of the target fibres: hyperpolarisation <= 0.3 h
EXCLUDE = [2]                          # the proximal stimulating contact
DESIGN_VARIANTS = {"free": False}     # name -> zero net current? "free" lets the current go anywhere (it
                                       # includes the balanced case, so it is never worse); add
                                       # "balanced": True to also see the best zero-net-current designs
TIME_LIMIT_s = 30                      # per mixed-integer solve (12-30 s each; ~6-7 min for the whole design)

T_TH_ON_s, T_PROX_ON_s = P["pulse"]["stimOnset_s"], 1.5e-3
PROX_WIDTH_s, T_AFTER_s = P["pulse"]["stimWidth_s"], P["pulse"]["tAfter_s"]
PLATEAU_START_s, PLATEAU_END_s = 0.1e-3, 0.9e-3
RISE_TIMES_s = [0.0]                   # rectangular (ramps never helped in rounds 2-3); add 0.3e-3 to test
RAMP_STEP_s = 10e-6
AMPS_A = np.geomspace(5e-6, 1e-3, 14)  # current of the largest-weight contact
CHECK_K = [2, 3]                       # designs checked with NEURON: every design with these K ...
EXTRA_CHECK = ["free_K1_N3"]           # ... plus these (a strong single-node design as a reference)
N_CHECK = 3                            # settings re-checked on the whole nerve
BLOCK_LIMIT = SIDE_LIMIT = RECRUIT
SHANNON_K = 1.75
N_WORKERS = max(1, (os.cpu_count() or 2) - 1)
''')

md(r'''
## 4. Data: unit activating functions, target fibres, maximum selective current

`anode_design.unit_af` interpolates the FEM potential of every contact (1 A) at the nodes of Ranvier of
every fibre line and takes the second difference: the matrix A (contacts × fibres × nodes). The thresholds
of the proximal pattern come from the cache of `block_volley.ipynb` (recomputed if the profiles changed);
the target fibres are the fibres of the target fascicle that fire at the maximum selective current.
''')

code(r'''
H5_PATH = BLOCK_DIR / f"{PATIENT}_fibre_profiles.h5"
Z_PROX_m = None
with h5py.File(H5_PATH, "r") as f:
    cz = dict(zip(f["contact_id"][:].ravel().astype(int).tolist(), f["contact_z_m"][:].ravel().tolist()))
Z_PROX_m = float(np.mean([cz[c] for c in PROX_PATTERN]))
ZN = anode_design.node_positions(FIBRE_DIAM_um, MRG["fibreLength_m"], Z_PROX_m)
t0 = time.time()
U = anode_design.unit_af(H5_PATH, ZN)
A, CID, fascicle = U["A"], U["contact_id"], U["fascicle"]
n_fibres = A.shape[1]
print(f"A: {A.shape} (contacts x fibres x nodes)  [{time.time()-t0:.0f} s]")

runner = fibre_jobs.FibreRunner(H5_PATH, MRG, FIBRE_DIAM_um, Z_PROX_m, MECH_DIR, N_WORKERS)
FIBRE = runner.fibre
NODE_B, NODE_M = runner.node_brain, runner.node_muscle
assert np.allclose(FIBRE.z_node_m, ZN), "node positions of the design and of the MRG fibre differ"

TH_CSV = BLOCK_DIR / "volley" / f"{PATIENT}_AS2_D{FIBRE_DIAM_um:g}_thresholds.csv"
TH_KEY = TH_CSV.with_name(TH_CSV.stem + "_key.json")
ok = TH_CSV.exists() and TH_KEY.exists() and json.loads(TH_KEY.read_text()).get("h5_mtime") == H5_PATH.stat().st_mtime
if ok and PROX_PATTERN == {2: -1.0}:
    th_A = pd.read_csv(TH_CSV)["threshold_A"].to_numpy()
    print("thresholds from the cache of block_volley.ipynb")
else:
    t_stop = T_TH_ON_s + PROX_WIDTH_s + T_AFTER_s
    jobs = [(i, PROX_PATTERN, T_TH_ON_s, PROX_WIDTH_s, t_stop, P["searchStart_A"], P["searchMax_A"], P["searchTol"])
            for i in range(n_fibres)]
    th_A = np.array(runner.map(fibre_jobs.threshold_job, jobs, chunksize=2))


def th10(k):
    t = np.sort(th_A[fascicle == k]); need = max(1, math.ceil(RECRUIT * len(t))); t = t[np.isfinite(t)]
    return t[need - 1] if len(t) >= need else np.nan


TH10 = {int(k): th10(k) for k in np.unique(fascicle)}
TARGET = TARGET_FASCICLE or int(min(TH10, key=lambda k: np.inf if np.isnan(TH10[k]) else TH10[k]))
I_PROX = SEL_MARGIN * min(v for k, v in TH10.items() if k != TARGET and np.isfinite(v))
TGT_ALL = np.flatnonzero(fascicle == TARGET)
TGT = np.flatnonzero((fascicle == TARGET) & (th_A <= I_PROX))
print(f"target F{TARGET}; I_prox = {I_PROX*1e6:.1f} uA fires {len(TGT)} of its {len(TGT_ALL)} fibres")
''')

md(r'''
## 5. Design

### 5.1 Where to block, and how much a perfect electrode could do

First with no limit on the number of contacts (a plain linear program): for each K and each possible
position of the K block nodes, the best ratio h. This is the upper bound of what this electrode can do.

**How to read this figure.** One curve per K: the best ratio h (block depth / worst depolarisation
anywhere) against the position of the first block node. Higher is better; the ratio falls as K grows
(more nodes to hyperpolarise at once). The single mirror anode AS 44 reached about 5 at K = 1 and never
blocked: so a design needs a clearly higher ratio, or more nodes at a similar ratio.
''')

code(r'''
allowed = [int(c) for c in CID if int(c) not in EXCLUDE]
zone = np.flatnonzero((ZN > BLOCK_ZONE_m[0]) & (ZN < BLOCK_ZONE_m[1]))
protect = np.flatnonzero(ZN < Z_PROX_m - 1e-3)
SCAN = []
t0 = time.time()
for vname, bal in DESIGN_VARIANTS.items():
    for K in K_LIST:
        for b in zone[:len(zone) - K + 1]:
            r = anode_design.design(A, CID, allowed=allowed, target_fibres=TGT, block_nodes=range(b, b + K),
                                    protect_nodes=protect, protect_frac=PROTECT_FRAC, balanced=bal)
            if r:
                SCAN.append({"variant": vname, "K": K, "b": int(b), "z_mm": ZN[b] * 1e3, "ratio": r["ratio"],
                             "n_active": r["n_active"]})
SCAN = pd.DataFrame(SCAN)
SCAN.to_csv(OUT_DIR / "scan_block_position.csv", index=False)
print(f"{len(SCAN)} linear programs in {time.time()-t0:.0f} s")
fig, axs = plt.subplots(1, len(DESIGN_VARIANTS), figsize=(13, 4), sharey=True, squeeze=False)
for ax, (vname, g) in zip(axs[0], SCAN.groupby("variant", sort=False)):
    for K, gk in g.groupby("K"):
        ax.plot(gk.z_mm, gk.ratio, ".-", label=f"K = {K}")
    ax.set(title=f"{vname}: all contacts allowed", xlabel="first block node z [mm]")
    ax.legend(fontsize=8)
axs[0, 0].set_ylabel("ratio h = block depth / worst depolarisation")
fig.savefig(OUT_DIR / "scan_block_position.png", dpi=130, bbox_inches="tight")
plt.show()
BEST_POS = SCAN.loc[SCAN.groupby(["variant", "K"]).ratio.idxmax()].reset_index(drop=True)
BEST_POS
''')

md(r'''
### 5.2 The fewest contacts: ratio against N

At the best position of each K, the mixed-integer program with at most N active contacts.

**How to read this figure.** For each variant and K, the best ratio reachable with at most N contacts.
The curve rises steeply for the first contacts and then flattens: the knee is the "minimum electrode".
The table lists every design (contacts and weights; anode > 0, cathode < 0, largest |weight| = 1).
''')

code(r'''
DESIGNS = []
t0 = time.time()
for _, row in BEST_POS.iterrows():
    bal = DESIGN_VARIANTS[row.variant]
    for N in N_LIST:
        r = anode_design.design(A, CID, allowed=allowed, target_fibres=TGT,
                                block_nodes=range(int(row.b), int(row.b) + int(row.K)), protect_nodes=protect,
                                protect_frac=PROTECT_FRAC, n_max=N, balanced=bal, time_limit_s=TIME_LIMIT_s)
        if r:
            name = f"{row.variant}_K{int(row.K)}_N{N}"
            DESIGNS.append({"name": name, "variant": row.variant, "K": int(row.K), "N_max": N,
                            "n_active": r["n_active"], "ratio": r["ratio"], "z_block_mm": float(row.z_mm),
                            "block_nodes": r["block_nodes"], "weights": r["weights"], "optimal": r["optimal"]})
DES = pd.DataFrame(DESIGNS)
print(f"{len(DES)} mixed-integer programs in {time.time()-t0:.0f} s")
fig, axs = plt.subplots(1, len(DESIGN_VARIANTS), figsize=(13, 4), sharey=True, squeeze=False)
for ax, (vname, g) in zip(axs[0], DES.groupby("variant", sort=False)):
    for K, gk in g.groupby("K"):
        ax.plot(gk.N_max, gk.ratio, "o-", label=f"K = {K}")
    ax.set(title=vname, xlabel="maximum number of active contacts N")
    ax.legend(fontsize=8)
axs[0, 0].set_ylabel("ratio h")
fig.savefig(OUT_DIR / "ratio_vs_contacts.png", dpi=130, bbox_inches="tight")
plt.show()
pd.set_option("display.max_colwidth", 120)
DES.assign(weights=DES.weights.map(lambda w: {c: round(v, 2) for c, v in w.items()}))[
    ["name", "n_active", "ratio", "z_block_mm", "weights"]]
''')

md(r'''
### 5.3 Candidates for the NEURON check

The designs to check: every design with K in `CHECK_K` (2 or 3 block nodes: one node is not enough, the
AP jumps over it) for every number of contacts, plus `EXTRA_CHECK` as a reference. Ratios of different K
are not comparable (more block nodes = lower ratio), so the choice is left to the NEURON check. Saved to `candidates.json` (also the input of
`block_volley.ipynb` with `DIST_PATTERN = "design:<name>"`).

**How to read this figure.** The AF of each candidate along the target fibres (median, per A of the
largest contact): the negative trough over the K block nodes is the block, everything positive must stay
small. Grey dashed: the sleeve ends; dotted: the two sites.
''')

code(r'''
CANDS = (DES[DES.K.isin(CHECK_K) | DES.name.isin(EXTRA_CHECK)]
         .assign(wkey=lambda d: d.weights.map(lambda w: tuple(sorted((c, round(v, 3)) for c, v in w.items()))))
         .drop_duplicates("wkey").drop(columns="wkey")          # same design found for two N: keep one
         .sort_values(["K", "n_active"]).reset_index(drop=True))
print(f"{len(CANDS)} designs to check:", list(CANDS.name))
cand_json = {r["name"]: {k: r[k] for k in ["variant", "K", "n_active", "ratio", "z_block_mm", "block_nodes"]}
             | {"weights": {str(c): w for c, w in r["weights"].items()}} for _, r in CANDS.iterrows()}
(OUT_DIR / "candidates.json").write_text(json.dumps(cand_json, indent=1, default=float))
fig, axs = plt.subplots(len(CANDS), 1, figsize=(11, 2.2 * len(CANDS)), sharex=True, squeeze=False)
for ax, (_, r) in zip(axs[:, 0], CANDS.iterrows()):
    af = anode_design.af_of(A, CID, r["weights"])
    ax.plot(ZN * 1e3, np.median(af[TGT], 0), ".-", color="tab:red", label=f"F{TARGET} target fibres (median)")
    ax.plot(ZN * 1e3, np.max(af[fascicle != TARGET], 0), color="grey", lw=0.8, label="other fascicles (max)")
    for s in (-25, 25):
        ax.axvline(s, color="k", ls="--", lw=0.6)
    for zz in (Z_PROX_m, np.mean(ZN[r["block_nodes"]])):
        ax.axvline(zz * 1e3, color="grey", ls=":", lw=0.8)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_title(f"{r['name']}: {r['n_active']} contacts {dict((c, round(w, 2)) for c, w in r['weights'].items())}, "
                 f"ratio {r['ratio']:.1f}", fontsize=8)
    ax.set_ylabel("AF [V/m^2 per A]", fontsize=8)
axs[0, 0].legend(fontsize=7)
axs[-1, 0].set_xlabel("z [mm]")
fig.tight_layout()
fig.savefig(OUT_DIR / "candidates_af.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
## 6. NEURON check on the target fascicle

For every candidate, rectangular and ramped anodal pulses (plateau 0.1–0.9 ms after the proximal onset),
currents `AMPS_A` on the largest-weight contact (the others scale with their weights). Two runs per fibre
of the target fascicle: proximal + distal pattern, and the distal pattern alone. Per setting:

* `frac_muscle` (C1 < 10 % of the fascicle), `frac_brain` (C0: the sensation must stay), `block_ratio`
  (blocked / fired, Tovbis et al. 2026), `frac_fired_by_pattern` (the distal pattern alone fires them).

First the reference without any distal pattern (all fibres, needed for C0–C3).
''')

code(r'''
def rect(t_on, w, a):
    return [(t_on, t_on + w, a)]


def trapezoid(t_on, t_off, a, rise, fall, step=RAMP_STEP_s):
    p = []
    n = int(round(rise / step))
    if n:
        e = t_on - rise + np.arange(n + 1) * rise / n
        p += [(e[k], e[k + 1], a * (k + 0.5) / n) for k in range(n)]
    p.append((t_on, t_off, a))
    n = int(round(fall / step))
    if n:
        e = t_off + np.arange(n + 1) * fall / n
        p += [(e[k], e[k + 1], a * (1 - (k + 0.5) / n)) for k in range(n)]
    return p


def dist_pulses(amp, rise):
    return trapezoid(T_PROX_ON_s + PLATEAU_START_s, T_PROX_ON_s + PLATEAU_END_s, amp, rise, rise)


def t_stop_for(p):
    return max(off for _, off, _ in p) + T_AFTER_s


def safety(weights, pulses):
    """Charge density per phase on the contact with the largest |weight| (all contacts are 0.2 x 0.2 mm)."""
    q = sum((off - on) * a for on, off, a in pulses) * max(abs(w) for w in weights.values())
    ok, Q_uC, D = volley.shannon_ok(q, 0.2e-3 * 0.2e-3, SHANNON_K)
    return {"charge_nC": Q_uC * 1e3, "density_uC_cm2": D, "shannon_ok": ok}


prox = (PROX_PATTERN, rect(T_PROX_ON_s, PROX_WIDTH_s, I_PROX))
t0 = time.time()
AP_BASE = runner.run_ap_times([(i, [prox], T_PROX_ON_s + PROX_WIDTH_s + T_AFTER_s) for i in range(n_fibres)])
OUT_BASE = volley.fascicle_outcomes(AP_BASE, fascicle, NODE_B, NODE_M)
FIRED = (volley.counts_at(AP_BASE, NODE_B) + volley.counts_at(AP_BASE, NODE_M)) > 0
print(f"reference: {FIRED[TGT_ALL].sum()} fibres of F{TARGET} fire  [{time.time()-t0:.0f} s]")
''')

code(r'''
rows, AP1 = [], {}
jobs, keys = [], []
for _, c in CANDS.iterrows():
    w = c["weights"]
    for rise in RISE_TIMES_s:
        for amp in AMPS_A:
            pul = dist_pulses(amp, rise)
            ts = t_stop_for(pul)
            jobs += [(int(i), [(w, pul)], ts) for i in TGT_ALL] + [(int(i), [prox, (w, pul)], ts) for i in TGT_ALL]
            keys.append((c["name"], rise, amp))
print(f"{len(jobs)} simulations on {N_WORKERS} processes ...")
t0 = time.time()
out = runner.run_ap_times(jobs)
print(f"done in {time.time()-t0:.0f} s")
nT = len(TGT_ALL)
fired_t = FIRED[TGT_ALL]
for s, (name, rise, amp) in enumerate(keys):
    alone, comb = out[2 * s * nT:(2 * s + 1) * nT], out[(2 * s + 1) * nT:(2 * s + 2) * nT]
    AP1[(name, rise, amp)] = (alone, comb)
    nm, nb = volley.counts_at(comb, NODE_M), volley.counts_at(comb, NODE_B)
    w = CANDS.set_index("name").loc[name, "weights"]
    rows.append({"name": name, "rise_ms": rise * 1e3, "amp_uA": amp * 1e6,
                 "frac_muscle": float(np.mean(nm > 0)), "frac_brain": float(np.mean(nb > 0)),
                 "block_ratio": float(np.mean(nm[fired_t] == 0)),
                 "frac_fired_by_pattern": float(np.mean((volley.counts_at(alone, NODE_B) + volley.counts_at(alone, NODE_M)) > 0)),
                 "started_at_fibre_end": float(np.mean(volley.started_at_end(alone, FIBRE.n_nodes))),
                 **safety(w, dist_pulses(amp, rise))})
S1 = pd.DataFrame(rows)
S1["C0"] = S1.frac_brain >= 0.9 * OUT_BASE.frac_brain[TARGET]
S1["C1"] = S1.frac_muscle < BLOCK_LIMIT
S1.to_csv(OUT_DIR / "check_target_fascicle.csv", index=False)
print("settings with C0 and C1 on the target fascicle:", int((S1.C0 & S1.C1).sum()))
S1[S1.C0 & S1.C1].round(3)
''')

md(r'''
### How to read this figure: block curves of each design

One panel per design, one line per ramp. *Solid*: fraction of the target fascicle still reaching the muscle
(C1: below the dotted line). *Dashed*: fraction fired by the distal pattern alone. *Dash-dot*: fraction
reaching the brain (C0: must stay at the level without the distal pattern, grey line). A working design
shows a range of currents where the solid line is below 10 % while the dashed line is still at 0.
''')

code(r'''
fig, axs = plt.subplots(1, len(CANDS), figsize=(4.2 * len(CANDS), 3.6), sharey=True, squeeze=False)
for ax, (_, c) in zip(axs[0], CANDS.iterrows()):
    g = S1[S1.name == c["name"]]
    for rise, gr in g.groupby("rise_ms"):
        (ln,) = ax.semilogx(gr.amp_uA, gr.frac_muscle, "-o", ms=3, label=f"ramps {rise:.1f} ms")
        ax.semilogx(gr.amp_uA, gr.frac_fired_by_pattern, "--", color=ln.get_color())
        ax.semilogx(gr.amp_uA, gr.frac_brain, "-.", color=ln.get_color(), lw=0.8)
    ax.axhline(BLOCK_LIMIT, color="k", ls=":", lw=0.8)
    ax.axhline(OUT_BASE.frac_brain[TARGET], color="grey", lw=0.6)
    ax.set_title(f"{c['name']} ({c['n_active']} contacts)", fontsize=8)
    ax.set_xlabel("current of the largest contact [uA]")
axs[0, 0].set_ylabel(f"fraction of F{TARGET}")
axs[0, 0].legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT_DIR / "check_target_curves.png", dpi=130, bbox_inches="tight")
plt.show()
''')

md(r'''
## 7. Whole-nerve check of the best settings

The `N_CHECK` settings that meet C0 and C1 on the target fascicle, fewest contacts first, then lowest
charge; each run on every fibre (distal pattern alone, and with the proximal pulse). Pass = C0–C3.
''')

code(r'''
ok1 = S1[S1.C0 & S1.C1].merge(CANDS[["name", "n_active"]], on="name")
SEL = ok1.sort_values(["n_active", "charge_nC"]).head(N_CHECK)
if SEL.empty:
    print("no design meets C0 and C1 on the target fascicle: nothing to check on the whole nerve")
CHECK, OUTF = [], {}
for _, s in SEL.iterrows():
    w = CANDS.set_index("name").loc[s["name"], "weights"]
    pul = dist_pulses(s.amp_uA * 1e-6, s.rise_ms * 1e-3)
    ts = t_stop_for(pul)
    t0 = time.time()
    res = runner.run_ap_times([(i, [(w, pul)], ts) for i in range(n_fibres)]
                              + [(i, [prox, (w, pul)], ts) for i in range(n_fibres)])
    alone, comb = res[:n_fibres], res[n_fibres:]
    o = volley.fascicle_outcomes(comb, fascicle, NODE_B, NODE_M, baseline=AP_BASE)
    OUTF[(s["name"], s.rise_ms, s.amp_uA)] = (o, alone, comb)
    oth = o.drop(index=TARGET)
    row = {"name": s["name"], "n_active": int(s.n_active), "rise_ms": s.rise_ms, "amp_uA": s.amp_uA,
           "C0": bool(o.frac_brain[TARGET] >= 0.9 * OUT_BASE.frac_brain[TARGET]),
           "C1": bool(o.frac_muscle[TARGET] < BLOCK_LIMIT),
           "C2": bool((oth.frac_muscle < SIDE_LIMIT).all()), "C3": bool((oth.frac_extra_brain < SIDE_LIMIT).all()),
           "target_frac_muscle": o.frac_muscle[TARGET], "target_block_ratio": o.block_ratio[TARGET],
           "charge_nC": s.charge_nC, "density_uC_cm2": s.density_uC_cm2, "shannon_ok": s.shannon_ok}
    row["PASS"] = row["C0"] and row["C1"] and row["C2"] and row["C3"]
    CHECK.append(row)
    o.to_csv(OUT_DIR / f"outcomes_{s['name']}_r{s.rise_ms:.1f}_{s.amp_uA:.0f}uA.csv")
    print(f"{s['name']} ramps {s.rise_ms:.1f} ms {s.amp_uA:.0f} uA: PASS={row['PASS']}  [{time.time()-t0:.0f} s]")
CHECK = pd.DataFrame(CHECK)
CHECK.to_csv(OUT_DIR / "check_whole_nerve.csv", index=False)
CHECK
''')

md(r'''
## 8. Result

The minimum configuration: the passing setting with the fewest active contacts, then the lowest charge.
To see it in motion: in `block_volley.ipynb` set `DIST_PATTERN = "design:<name>"`, run it, and export it
to the viewer (its section 10).
''')

code(r'''
best = CHECK[CHECK.PASS].sort_values(["n_active", "charge_nC"]).head(1) if len(CHECK) else CHECK
summary = {"patient": PATIENT, "target": TARGET, "I_prox_uA": I_PROX * 1e6, "fired_target": int(len(TGT)),
           "designs": DES.to_dict(orient="records"), "candidates": list(CANDS.name),
           "check_target": S1.to_dict(orient="records"),
           "check_whole_nerve": CHECK.to_dict(orient="records") if len(CHECK) else [],
           "best": best.to_dict(orient="records")}
(OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
if len(best):
    b = best.iloc[0]
    w = CANDS.set_index("name").loc[b["name"], "weights"]
    print(f"MINIMUM CONFIGURATION: {b['name']} -> {b.n_active} contacts {w}, {b.amp_uA:.0f} uA on the largest, "
          f"ramps {b.rise_ms:.1f} ms; charge {b.charge_nC:.1f} nC, {b.density_uC_cm2:.0f} uC/cm^2 (Shannon ok: {b.shannon_ok})")
else:
    print("no configuration passes C0-C3: see the tables above for how close each design gets")
runner.close()
''')

nb = {"cells": [], "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
for i, (kind, src) in enumerate(cells):
    lines = src.split("\n")
    cell = {"cell_type": kind, "id": f"d{i:02d}", "metadata": {}, "source": [l + "\n" for l in lines[:-1]] + [lines[-1]]}
    if kind == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    nb["cells"].append(cell)
OUT.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print("wrote", OUT, len(cells), "cells")
