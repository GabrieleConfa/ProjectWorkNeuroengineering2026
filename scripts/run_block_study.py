"""Conduction-block study: can the distal site stop the APs the proximal site sends towards the hand?

What it does (docs/block.md for the walkthrough, docs/decisions.md D-12 for the choices)
    Reads the FEM potential along the fibres (results/block/<id>_fibre_profiles.h5, written by
    scripts/export_block_profiles.m) and the parameters (results/block/<id>_block_params.json, a copy
    of config/block_params.m), and runs the MRG fibre model (src/block/mrg_fibre.py) for one pattern
    (proximal contacts + weights + the fascicle it is meant to excite):

    0. checks: no AP without stimulus; conduction velocity of the fibre.
    1. thresholds of the TARGET fascicle: smallest cathodic current on the proximal pattern (and on
       the distal mirror pattern) that launches an AP, for each of its fibres; the 10%-of-fibres
       threshold (Zelechowski 2020) sets I_stim = stimFactor x that threshold. Every other fibre is
       simulated once at I_stim ("does it fire?"); --all-thresholds computes all of them.
    2. proximal alone / distal alone at I_stim, on all fibres: which fibres each site excites.
    3. collision: the distal site repeats the same pattern (same polarity and amplitude) after a delay.
    4. cathodal block (supervisor's scheme): the same distal electrodes as cathodes, with a larger
       amplitude, swept over amplitude x delay x pulse width on the target fascicle; the best settings
       are re-checked on all fibres (off-target excitation by the strong distal cathode).
    5. anodal block: the distal contacts with reversed polarity (anodic), 1 ms, and the anode alone.

    For every fibre, the expected outcome is: as many APs at the proximal end (brain) as the proximal
    pulse alone gives (1 if the fibre is excited, else 0) and NO AP at the distal end (muscle).
    APs are counted at the nodes detectInset_m from each fibre end. Proximal = z < 0.

Outputs (results/block/)
    <tag>_thresholds.csv     per fibre: position, fascicle, thresholds (target fascicle), fires at I_stim
    <tag>_scenarios.csv      per fibre and scenario: APs at the proximal/distal end
    <tag>_summary.json       the numbers printed at the end
    <tag>_raster.png         membrane voltage along one fibre vs time, per scenario
    <tag>_sweeps.png         collision vs delay, anodal block vs amplitude
    <tag>_cathodal.png       success of the cathodal block vs amplitude and delay, per pulse width
    (<tag> = <id>_<pattern>_D<diameter>)

Run (from the project root, after export_block_profiles in MATLAB):
    python scripts/run_block_study.py --pattern AS2 [--diam 10.0] [--workers 7] [--all-thresholds]
Time: ~0.3 s per fibre simulation, a few thousand simulations: a few minutes on 7 processes.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "block"))
from mrg_fibre import MRG_TABLE, Channel, MRGFibre  # noqa: E402

MECH_DIR = ROOT / "results" / "mrg" / "model"
OUT_DIR = ROOT / "results" / "block"

# ---------------------------------------------------------------------------------------------------
# Worker side: every process loads the profiles once and keeps one fibre object per z alignment.
_W: dict = {}


def _init_worker(h5_path: str, params: dict, diam_um: float):
    """Run once in every worker process: load the FEM profiles into memory and remember the parameters.

    h5_path  results/block/<id>_fibre_profiles.h5;  params  the block parameters (JSON copy of
    config/block_params.m);  diam_um  fibre diameter [um]. The MRG fibres are built lazily by _fibre.
    """
    with h5py.File(h5_path, "r") as f:                 # everything is copied into memory ([:]): the file is closed
        _W["z"] = f["z_m"][:].ravel()                  # sample points along the nerve [m]
        _W["V"] = f["V_perA"][:]                       # (nC, nL, nZ) float32 [V per A], ~25 MB
        _W["cid"] = list(f["contact_id"][:].ravel().astype(int))
    _W["p"], _W["diam"], _W["fibres"] = params, diam_um, {}


def _fibre(z_align_m: float) -> MRGFibre:
    """MRG fibre (NEURON) whose node of Ranvier sits at z_align_m [m]; built once per alignment and reused.

    With nodeAlignment = "proximalStation" all fibres share one alignment, so only one fibre is built.
    """
    key = round(z_align_m * 1e9)                       # alignment rounded to 1 nm as dictionary key
    if key not in _W["fibres"]:
        m = _W["p"]["mrg"]
        _W["fibres"][key] = MRGFibre(_W["diam"], m["fibreLength_m"], z_align_m, MECH_DIR, m["temperature_C"])
    return _W["fibres"][key]


def _profile(i: int, pattern: dict[int, float]) -> np.ndarray:
    """Potential along fibre i for 1 A of pattern current: sum of weight x contact profile."""
    return sum(w * _W["V"][_W["cid"].index(c), i, :].astype(float) for c, w in pattern.items())


def _simulate(i, z_align_m, chans, t_stop_s, record=False):
    """chans = list of (pattern dict, pulses list). Returns (APs prox end, APs dist end, result, fibre)."""
    m = _W["p"]["mrg"]
    fib = _fibre(z_align_m)
    channels = [Channel(_profile(i, pat), pulses) for pat, pulses in chans]
    res = fib.run(channels, _W["z"], t_stop_s, m["dt_s"], m["vInit_V"], m["apThreshold_V"], record)
    half = m["fibreLength_m"] / 2
    n_p = len(res["ap_times_s"][fib.node_index(-half + m["detectInset_m"])])
    n_d = len(res["ap_times_s"][fib.node_index(half - m["detectInset_m"])])
    return n_p, n_d, res, fib


def task_threshold(args):
    """Smallest cathodic amplitude of `pattern` that sends an AP to either end (NaN if > max)."""
    i, z_align_m, pattern, width_s, onset_s, t_stop_s = args
    p = _W["p"]
    fires = lambda I: sum(_simulate(i, z_align_m, [(pattern, [(onset_s, onset_s + width_s, I)])],
                                    t_stop_s)[:2]) > 0
    lo, hi = 0.0, p["searchStart_A"]
    while not fires(hi):
        lo, hi = hi, hi * 2
        if lo >= p["searchMax_A"]:
            return i, float("nan")
    while (hi - lo) > p["searchTol"] * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if fires(mid) else (mid, hi)
    return i, hi


def task_scenario(args):
    """One fibre under one scenario: returns (i, label, APs at proximal end, APs at distal end)."""
    i, z_align_m, label, chans, t_stop_s = args
    n_p, n_d, _, _ = _simulate(i, z_align_m, chans, t_stop_s)
    return i, label, n_p, n_d


def task_record(args):
    """One fibre with the membrane voltage recorded at every node (for the raster plots).

    args = (i, z_align_m, chans, t_stop_s); returns (t_s, v_V (nodes x time), z of the nodes [m]).
    """
    i, z_align_m, chans, t_stop_s = args
    _, _, res, fib = _simulate(i, z_align_m, chans, t_stop_s, record=True)
    return res["t_s"], res["v_V"], fib.z_node_m


# ---------------------------------------------------------------------------------------------------
def main():
    """Parse the command line, run steps 0-5 of the module docstring on a process pool, write the outputs."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patient", default="UA01")
    ap.add_argument("--pattern", default=None, help="name of a pattern in block_params.m (default: first)")
    ap.add_argument("--diam", type=float, default=None, help="fibre diameter [um] (default: b.mrg.fibreDiam_um)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--all-thresholds", action="store_true",
                    help="exact thresholds for every fibre (slow: ~15 min), not only the target fascicle")
    args = ap.parse_args()

    h5_path = OUT_DIR / f"{args.patient}_fibre_profiles.h5"
    p = json.loads((OUT_DIR / f"{args.patient}_block_params.json").read_text())
    m, pu, cb = p["mrg"], p["pulse"], p["cathodal"]
    diam = args.diam or m["fibreDiam_um"]
    pats = p["patterns"] if isinstance(p["patterns"], list) else [p["patterns"]]
    pat = next(q for q in pats if q["name"] == args.pattern) if args.pattern else pats[0]
    nw = args.workers or p["nWorkers"] or max(1, (os.cpu_count() or 2) - 1)

    with h5py.File(h5_path, "r") as f:
        z = f["z_m"][:].ravel()
        x, y = f["x_m"][:].ravel(), f["y_m"][:].ravel()
        fas = f["fascicle"][:].ravel().astype(int)
        cid = f["contact_id"][:].ravel().astype(int)
        mfrom = f["mirrored_from"][:].ravel().astype(int)
        cz = f["contact_z_m"][:].ravel()
        source = f.attrs["source"]
    nL, nF = len(x), int(fas.max())

    prox = {int(c): float(w) for c, w in zip(np.atleast_1d(pat["contacts"]), np.atleast_1d(pat["weights"]))}
    partner = {int(s): int(c) for c, s in zip(cid, mfrom) if s > 0}
    partner.update({v: k for k, v in partner.items()})
    dist = {partner[c]: w for c, w in prox.items()}      # same electrodes on the other site, same polarity
    anod = {c: -w for c, w in dist.items()}              # same electrodes, reversed polarity
    zc = dict(zip(cid.tolist(), cz.tolist()))
    assert all(zc[c] < 0 for c in prox) and all(zc[c] > 0 for c in dist), "proximal = z < 0 (D-12)"
    z_prox, z_dist = zc[next(iter(prox))], zc[next(iter(dist))]
    gap = z_dist - z_prox
    target = int(pat["target"])
    tgt = np.flatnonzero(fas == target)

    if m["nodeAlignment"] == "random":
        rng = np.random.default_rng(m["seed"])
        zal = z_prox + rng.uniform(0, MRG_TABLE[diam][5] * 1e-6, nL)
    else:
        zal = np.full(nL, z_prox)

    tag = f"{args.patient}_{pat['name']}_D{diam:g}"
    print(f"Pattern {pat['name']}: proximal {prox} at z = {z_prox*1e3:+.1f} mm, distal {dist} at "
          f"z = {z_dist*1e3:+.1f} mm (gap {gap*1e3:.1f} mm), target F{target} ({len(tgt)} fibres). "
          f"MRG D = {diam} um, {nL} fibres, {nw} processes.")
    print(f"FEM source: {source}")
    on, w_s, t_after = pu["stimOnset_s"], pu["stimWidth_s"], pu["tAfter_s"]
    t_single = on + w_s + t_after
    t0 = time.time()
    out = {}                                             # (label, fibre) -> (APs brain, APs muscle)

    def run_scenarios(ex, scen, fibres):
        """Run every scenario (label, channels, t_stop) on every fibre in parallel; results go to `out`."""
        jobs =[(int(i), zal[i], lab, ch, ts) for lab, ch, ts in scen for i in fibres]
        for i, lab, n_p, n_d in ex.map(task_scenario, jobs, chunksize=4):
            out[(lab, i)] = (n_p, n_d)

    with ProcessPoolExecutor(nw, initializer=_init_worker, initargs=(str(h5_path), p, diam)) as ex:
        # 0. checks ----------------------------------------------------------------------------------
        t, v, zn = ex.submit(task_record, (0, zal[0], [], t_single)).result()
        assert np.all(v < m["apThreshold_V"]), "spontaneous AP without stimulus"

        # 1. thresholds (target fascicle, or all with --all-thresholds), cached by their inputs --------
        th_csv, th_key_file = OUT_DIR / f"{tag}_thresholds.csv", OUT_DIR / f"{tag}_thresholds_key.json"
        th_fibres = np.arange(nL) if args.all_thresholds else tgt
        th_key = {"source": str(source), "nL": nL, "z": [float(z[0]), float(z[-1]), len(z)],
                  "fibres": [p["gridStepXY_m"], p["fibreStepXY_m"], p["fibreInset_m"]],
                  "mrg": m, "pulse": [w_s, on, t_after], "contacts": [prox, dist],
                  "search": [p["searchStart_A"], p["searchMax_A"], p["searchTol"]],
                  "fibresWithThreshold": th_fibres.tolist()}
        th_key = json.loads(json.dumps(th_key))
        cached = th_key_file.exists() and th_csv.exists() and \
            json.loads(th_key_file.read_text()) == th_key
        if cached:
            th = np.genfromtxt(th_csv, delimiter=",", names=True)
            thP, thD = th["thProx_uA"] * 1e-6, th["thDist_uA"] * 1e-6
            print(f"1. thresholds reused from {th_csv.name} (same inputs)")
        else:
            thP, thD = np.full(nL, np.nan), np.full(nL, np.nan)
            for pattern, arr in [(prox, thP), (dist, thD)]:
                jobs = [(int(i), zal[i], pattern, w_s, on, t_single) for i in th_fibres]
                for i, th in ex.map(task_threshold, jobs, chunksize=2):
                    arr[i] = th
            print(f"1. thresholds of {len(th_fibres)} fibres, proximal and distal [{time.time()-t0:.0f} s]")

        # 10%-of-fibres threshold of every fascicle whose fibres all have a computed threshold
        # (NaN threshold of a computed fibre = no AP up to searchMax_A)
        computed = np.isin(np.arange(nL), th_fibres)
        th10 = np.full(nF + 1, np.nan)
        for k in range(1, nF + 1):
            if not np.all(computed[fas == k]):
                continue
            a = np.sort(thP[fas == k][np.isfinite(thP[fas == k])])
            n_need = max(1, int(np.ceil(p["recruitFraction"] * np.sum(fas == k))))
            if len(a) >= n_need:
                th10[k] = a[n_need - 1]
        assert np.isfinite(th10[target]), f"target F{target} not excited up to {p['searchMax_A']*1e6:.0f} uA"
        I_stim = p["stimFactor"] * th10[target]
        print(f"   F{target} (target): 10%-of-fibres threshold {th10[target]*1e6:.1f} uA = "
              f"{th10[target]*w_s*1e9:.2f} nC ({w_s*1e6:.0f} us); I_stim = {p['stimFactor']} x = "
              f"{I_stim*1e6:.1f} uA ({I_stim*w_s*1e9:.2f} nC per phase)")
        if args.all_thresholds:
            for k in range(1, nF + 1):
                if k != target:
                    print(f"     F{k:<2d} {th10[k]*1e6:7.1f} uA = {th10[k]*w_s*1e9:5.2f} nC")
        ok = np.isfinite(thP) & np.isfinite(thD)
        if ok.any():
            r = thD[ok] / thP[ok]
            print(f"   mirror check, distal / proximal threshold of the same fibre: median {np.median(r):.2f}, "
                  f"range {r.min():.2f} ... {r.max():.2f}")

        # 2. which fibres each site excites at I_stim (all fibres, one run each) -----------------------
        stim = (prox, [(on, on + w_s, I_stim)])
        run_scenarios(ex, [("proximal only", [stim], t_single),
                           ("distal only", [(dist, [(on, on + w_s, I_stim)])], t_single)], range(nL))
        recP = np.array([sum(out[("proximal only", i)]) > 0 for i in range(nL)])
        recD = np.array([sum(out[("distal only", i)]) > 0 for i in range(nL)])
        expect_brain = np.array([out[("proximal only", i)][0] for i in range(nL)])
        fib_set = np.flatnonzero(recP | recD)
        print(f"2. at I_stim: proximal excites {recP.sum()} fibres "
              f"({', '.join(f'F{k}: {np.sum(recP & (fas == k))}' for k in range(1, nF + 1) if np.any(recP & (fas == k)))}), "
              f"distal {recD.sum()} [{time.time()-t0:.0f} s]")

        # 3. collision: same pattern, same amplitude, delayed --------------------------------------------
        coll = [(f"collision delay {d*1e3:.2f} ms", [stim, (dist, [(on + d, on + d + w_s, I_stim)])],
                 on + d + w_s + t_after) for d in pu["delaySweep_s"]]
        run_scenarios(ex, coll, fib_set)

        # 4. cathodal block: same electrodes as cathodes, larger amplitude -------------------------------
        cset = np.union1d(tgt, fib_set)
        cath = []
        for wd in cb["widths_s"]:
            for d in cb["delays_s"]:
                for fa in cb["factors"]:
                    lab = f"cathodal x{fa:g} w{wd*1e6:.0f}us d{d*1e3:.1f}ms"
                    cath.append((lab, [stim, (dist, [(on + d, on + d + wd, fa * I_stim)])],
                                 on + d + wd + t_after, fa, wd, d))
        run_scenarios(ex, [c[:3] for c in cath], cset)
        print(f"3-4. collision ({len(coll)} delays) on {len(fib_set)} fibres, cathodal grid ({len(cath)}) on "
              f"{len(cset)} fibres [{time.time()-t0:.0f} s]")

        # 5. anodal block: reversed polarity, and the anode alone ---------------------------------------
        wa = pu["anodalWidth_s"]
        anodal = [(f"anodal x{fa:g}", [stim, (anod, [(on, on + wa, fa * I_stim)])], on + wa + t_after)
                  for fa in p["anodalFactors"]]
        anodal += [(f"anode alone x{fa:g}", [(anod, [(on, on + wa, fa * I_stim)])], on + wa + t_after)
                   for fa in p["anodalFactors"]]
        run_scenarios(ex, anodal, fib_set)

        # ---- success of each scenario ----------------------------------------------------------------
        def stats(lab, fibres):
            """Percentages of `fibres` for scenario `lab`: AP to the muscle, brain APs as expected / extra /
            lost, success (no AP to the muscle AND the brain receives what the proximal pulse alone gives)."""
            fibres =[int(i) for i in fibres if (lab, int(i)) in out]
            r = np.array([out[(lab, i)] for i in fibres]).reshape(-1, 2)
            e = expect_brain[fibres]
            n = len(fibres)
            pct = lambda b: float(100 * np.mean(b)) if n else float("nan")
            return {"n": n, "muscle_pct": pct(r[:, 1] > 0), "brain_ok_pct": pct(r[:, 0] == e),
                    "brain_extra_pct": pct(r[:, 0] > e), "brain_lost_pct": pct(r[:, 0] < e),
                    "success_pct": pct((r[:, 1] == 0) & (r[:, 0] == e))}

        # best cathodal settings on the target set, then re-checked on all fibres
        cstat = [(stats(c[0], cset), c) for c in cath]
        best = sorted(cstat, key=lambda sc: (-sc[0]["success_pct"], sc[1][3] * sc[1][4]))[:cb["nCheckAll"]]
        best = [b for b in best if b[0]["success_pct"] > 0]
        if best:
            run_scenarios(ex, [b[1][:3] for b in best], range(nL))
        print(f"   cathodal: best success on the target set {best[0][0]['success_pct']:.0f}%" if best else
              "   cathodal: no setting succeeds on any fibre")

        # example fibre and rasters
        cand = np.flatnonzero((fas == target) & recP)
        iex = int(cand[np.argmin(np.where(np.isfinite(thP[cand]), thP[cand], np.inf))])
        # collision delay shown in the raster: the swept delay closest to half the travel time over the gap,
        # with the rule of thumb CV ~ 5 m/s per um of diameter (about 50 m/s for the 10 um MRG fibre)
        d_ex = min(pu["delaySweep_s"], key=lambda d: abs(d - 0.5 * gap / (5.0 * diam)))
        rasters = [("proximal only", [stim], t_single),
                   (f"collision: distal repeat after {d_ex*1e3:.2f} ms",
                    [stim, (dist, [(on + d_ex, on + d_ex + w_s, I_stim)])], on + d_ex + w_s + t_after)]
        if best:
            _, (lab, ch, ts, fa, wd, d) = best[0]
            rasters.append((f"cathodal block: distal x{fa:g} ({fa*I_stim*1e6:.0f} uA), "
                            f"{wd*1e6:.0f} us, after {d*1e3:.1f} ms", ch, ts))
        else:   # show the strongest, longest cathodal setting
            lab, ch, ts, fa, wd, d = max(cath, key=lambda c: (c[3], c[4]))
            rasters.append((f"cathodal x{fa:g} ({fa*I_stim*1e6:.0f} uA), {wd*1e6:.0f} us, "
                            f"after {d*1e3:.1f} ms (no success)", ch, ts))
        rasters.append((f"anodal: distal reversed x16 ({16*I_stim*1e6:.0f} uA, {wa*1e3:.1f} ms)",
                        [stim, (anod, [(on, on + wa, 16 * I_stim)])], on + wa + t_after))
        rec = [ex.submit(task_record, (iex, zal[iex], ch, ts)).result() for _, ch, ts in rasters]

    # conduction velocity from the proximal-only raster (distal half of the fibre)
    t_r, v_r, zn = rec[0]
    cross = lambda k: t_r[1:][(v_r[k, :-1] < 0) & (v_r[k, 1:] >= 0)]
    a, b_ = np.argmin(abs(zn - 5e-3)), np.argmin(abs(zn - 20e-3))
    cv = (zn[b_] - zn[a]) / (cross(b_)[0] - cross(a)[0]) if len(cross(a)) and len(cross(b_)) else np.nan
    print(f"0. checks: no AP without stimulus; conduction velocity {cv:.1f} m/s (D = {diam} um); "
          f"travel time over the {gap*1e3:.1f} mm gap {gap/cv*1e3:.2f} ms")

    # ---- report ------------------------------------------------------------------------------------
    summary = {"pattern": pat["name"], "contacts": [prox, dist], "target": target, "diam_um": diam,
               "I_stim_uA": I_stim * 1e6, "th10_target_uA": th10[target] * 1e6, "cv_mps": cv,
               "gap_mm": gap * 1e3, "n_excited_prox": int(recP.sum()), "n_excited_dist": int(recD.sum()),
               "scenarios": {}, "cathodal_best": []}
    hdr = f"{'scenario':40s} {'fibres':>6s} {'success':>8s} {'muscle':>7s} {'brain ok':>9s} {'extra':>6s} {'lost':>5s}"
    print(f"\nPer scenario, % of the fibres: success = no AP to the muscle and the same APs to the brain as "
          f"the proximal pulse alone\n{hdr}")

    def line(lab, s):
        """One row of the printed table (output of stats)."""
        print(f"{lab:40s} {s['n']:6d} {s['success_pct']:7.0f}% {s['muscle_pct']:6.0f}% {s['brain_ok_pct']:8.0f}% "
              f"{s['brain_extra_pct']:5.0f}% {s['brain_lost_pct']:4.0f}%")

    for lab, _, _ in [("proximal only", 0, 0), ("distal only", 0, 0)] + coll + anodal:
        s = stats(lab, fib_set)
        summary["scenarios"][lab] = s
        line(lab, s)
    print(f"cathodal grid on the target set ({len(cset)} fibres), best {len(best)}; then on all {nL} fibres:")
    for s_t, c in best:
        s_all = stats(c[0], range(nL))
        summary["cathodal_best"].append({"label": c[0], "factor": c[3], "width_us": c[4] * 1e6,
                                         "delay_ms": c[5] * 1e3, "I_uA": c[3] * I_stim * 1e6,
                                         "charge_nC": c[3] * I_stim * c[4] * 1e9, "target_set": s_t,
                                         "all_fibres": s_all})
        line(c[0] + " (target)", s_t)
        line(c[0] + " (all)", s_all)
    summary["cathodal_grid"] = {c[0]: s for s, c in cstat}

    with open(OUT_DIR / f"{tag}_thresholds.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["fibre", "x_mm", "y_mm", "fascicle", "thProx_uA", "thDist_uA", "firesProx", "firesDist"])
        for i in range(nL):
            wr.writerow([i, f"{x[i]*1e3:.3f}", f"{y[i]*1e3:.3f}", fas[i], f"{thP[i]*1e6:.4f}",
                         f"{thD[i]*1e6:.4f}", int(recP[i]), int(recD[i])])
    th_key_file.write_text(json.dumps(th_key, indent=1))
    with open(OUT_DIR / f"{tag}_scenarios.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["scenario", "fibre", "fascicle", "APs_proximal_end", "APs_distal_end"])
        for (lab, i), (n_p, n_d) in out.items():
            wr.writerow([lab, i, fas[i], n_p, n_d])
    (OUT_DIR / f"{tag}_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    plot_raster(OUT_DIR / f"{tag}_raster.png", rec, [r[0] for r in rasters], z_prox, z_dist, iex, fas[iex], diam)
    plot_sweeps(OUT_DIR / f"{tag}_sweeps.png", summary, pu, p, len(fib_set), gap, cv)
    plot_cathodal(OUT_DIR / f"{tag}_cathodal.png", cstat, cb, I_stim, len(cset), target, gap, cv)
    print(f"\nWrote {OUT_DIR / tag}_*.csv/json/png  [{time.time()-t0:.0f} s]")


# ---------------------------------------------------------------------------------------------------
def _plt():
    """matplotlib.pyplot with the non-interactive Agg backend (figures are only saved to files)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_raster(path, rec, titles, z_prox, z_dist, i, fasc, diam):
    """<tag>_raster.png: membrane voltage (colour) vs time (x) and node position (y), one panel per scenario.

    rec  list of (t_s, v_V, z_node_m) from task_record;  titles  one per panel;  z_prox, z_dist  the two
    sites [m] (dashed lines);  i, fasc  fibre number and its fascicle;  diam  fibre diameter [um].
    """
    plt = _plt()
    fig, axs = plt.subplots(1, len(rec), figsize=(4.6 * len(rec), 4.6), sharey=True)
    for ax, (t, v, zn), ttl in zip(np.atleast_1d(axs), rec, titles):
        im = ax.pcolormesh(t * 1e3, zn * 1e3, v * 1e3, cmap="magma", vmin=-100, vmax=40, shading="auto")
        for zz, lab in [(z_prox, "proximal site"), (z_dist, "distal site")]:
            ax.axhline(zz * 1e3, color="c", lw=0.8, ls="--")
            ax.text(t[-1] * 1e3, zz * 1e3, lab, color="c", ha="right", va="bottom", fontsize=8)
        ax.set_title(ttl, fontsize=8.5)
        ax.set_xlabel("time [ms]")
    np.atleast_1d(axs)[0].set_ylabel("node position along the nerve [mm]  (brain = bottom, muscle = top)")
    fig.colorbar(im, ax=axs, label="membrane potential at the nodes [mV]")
    fig.suptitle(f"MRG {diam} um fibre no. {i} (F{fasc}): APs are the bright diagonal stripes "
                 "(slope = conduction velocity)", fontsize=10)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_sweeps(path, summary, pu, p, n, gap, cv):
    """<tag>_sweeps.png: left, collision outcome vs delay of the distal repeat; right, anodal block vs current.

    summary  the summary dict (scenario statistics);  pu, p  pulse block and all parameters;  n  number of
    excited fibres;  gap  distance between the sites [m];  cv  conduction velocity [m/s] (dotted line).
    """
    plt = _plt()
    S = summary["scenarios"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    # left: collision scenarios, one point per delay of the distal repeat
    d = np.array(pu["delaySweep_s"]) * 1e3
    lab = [f"collision delay {x:.2f} ms" for x in d]
    a1.plot(d, [S[l]["muscle_pct"] for l in lab], "o-", label="AP reaches the muscle (want 0)")
    a1.plot(d, [S[l]["brain_ok_pct"] for l in lab], "s-", label="brain gets the right APs (want 100)")
    a1.plot(d, [S[l]["brain_extra_pct"] for l in lab], "^-", label="extra APs to the brain")
    if np.isfinite(cv):
        a1.axvline(gap / cv * 1e3, color="k", ls=":", lw=1)
        a1.text(gap / cv * 1e3, 30, " AP arrives at\n the distal site", fontsize=8)
    a1.set(xlabel="delay of the distal repeat [ms]", ylabel=f"% of the {n} excited fibres",
           title="Collision: same pattern repeated distally", ylim=(-5, 105))
    a1.legend(fontsize=8, loc="center right")
    # right: anodal block scenarios against the anodal current (factor x I_stim), and the anode alone
    fa = np.array(p["anodalFactors"], dtype=float)
    lab = [f"anodal x{x:g}" for x in fa]
    a2.semilogx(fa * summary["I_stim_uA"], [S[l]["muscle_pct"] for l in lab], "o-", label="AP reaches the muscle")
    a2.semilogx(fa * summary["I_stim_uA"], [S[l]["brain_ok_pct"] for l in lab], "s-", label="brain gets the right APs")
    a2.semilogx(fa * summary["I_stim_uA"], [S[f"anode alone x{x:g}"]["muscle_pct"] for x in fa], "k--",
                label="anode ALONE excites (AP to the muscle)")
    a2.set(xlabel="distal anodal current [uA]", ylabel="%", ylim=(-5, 105),
           title=f"Anodal block: distal contacts reversed, {pu['anodalWidth_s']*1e3:.1f} ms pulse")
    a2.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_cathodal(path, cstat, cb, I_stim, n, target, gap, cv):
    """<tag>_cathodal.png: success of the cathodal block, current (y) x onset delay (x), one panel per width.

    cstat  list of (stats, scenario tuple);  cb  the "cathodal" parameter block;  I_stim  proximal current
    [A];  n  fibres in the target set;  target  fascicle number;  gap [m], cv [m/s] for the arrival line.
    """
    plt = _plt()
    W, D, Fa = cb["widths_s"], cb["delays_s"], cb["factors"]
    grid = {(c[4], c[5], c[3]): s["success_pct"] for s, c in cstat}
    fig, axs = plt.subplots(1, len(W), figsize=(4.4 * len(W), 3.9), sharey=True)
    for ax, wd in zip(np.atleast_1d(axs), W):
        M = np.array([[grid[(wd, d, fa)] for d in D] for fa in Fa])
        im = ax.imshow(M, origin="lower", aspect="auto", vmin=0, vmax=100, cmap="viridis",
                       extent=[D[0] * 1e3 - 0.05, D[-1] * 1e3 + 0.05, -0.5, len(Fa) - 0.5])
        ax.set_yticks(range(len(Fa)), [f"{fa * I_stim * 1e6:.0f}" for fa in Fa])
        if np.isfinite(cv):
            ax.axvline(gap / cv * 1e3, color="w", ls=":", lw=1)
        ax.set(xlabel="distal onset after the proximal one [ms]",
               title=f"distal cathode, {wd * 1e6:.0f} us pulse")
    np.atleast_1d(axs)[0].set_ylabel("distal cathodic current [uA]")
    fig.colorbar(im, ax=axs, label=f"% success on the {n} fibres of F{target} / excited")
    fig.suptitle("Cathodal block (same electrodes as cathodes on the distal site): success = no AP to the "
                 "muscle, sensation unchanged; dotted = AP arrival", fontsize=9.5)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
