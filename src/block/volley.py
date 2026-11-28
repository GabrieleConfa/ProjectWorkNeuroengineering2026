"""Population ("volley") analysis of the conduction-block study.

Why this module exists
    The stimulator acts on, and can only observe, whole fascicles (D-15). The MRG simulations, however,
    are per fibre. This module turns the per-fibre results into fascicle-level quantities:

      * the VOLLEY of a fascicle: at each time t and each node position z, the fraction of the
        fascicle's fibres that carry an action potential (AP) there. Plotted against z it is a bump that
        travels along the nerve: the compound activity of the fascicle (what a recording electrode
        around it would sum up);
      * the REACH profile: for each node, the fraction of the fascicle's fibres in which an AP got at
        least as far as that node. Without a block it is flat; an anodal block shows up as a drop at
        the distal site;
      * the OUTCOMES per fascicle (fibres that delivered an AP to the brain end and to the muscle end),
        and the block scores of docs/decisions.md D-16.

    It works only on the output of fibre_jobs.ap_times_job: for every simulated fibre, a list with one
    array per node of Ranvier holding the times [s] of the upward crossings of the AP threshold (0 mV).
    It does not run NEURON, so it can also be used to re-analyse saved results.

Conventions
    ap      list over fibres (index = fibre line in the profile file) of lists over nodes of arrays of
            crossing times [s]. Node 0 is the proximal end (towards the brain, z < 0).
    SI units in the code (m, s); mm / ms only in labels.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------------------------------
# Counting
# ----------------------------------------------------------------------------------------------------
def counts_at(ap: list, node: int) -> np.ndarray:
    """Number of APs recorded at `node`, for every fibre of `ap` (array of ints)."""
    return np.array([len(a[node]) for a in ap], dtype=int)


def fascicle_outcomes(ap: list, fascicle: np.ndarray, node_brain: int, node_muscle: int,
                      baseline: list | None = None) -> pd.DataFrame:
    """Per-fascicle table of what reaches the two ends of the fibres.

    ap          AP times of one scenario, for ALL fibres (len = number of fibre lines)
    fascicle    fascicle number of every fibre line (1 ... nF)
    node_brain, node_muscle   the detection nodes (5 mm from each fibre end)
    baseline    optional: AP times of the reference scenario (proximal pulse alone). When given, the
                table also says which fibres the proximal pulse fired and what happened to them.

    Columns (one row per fascicle, "frac_" = fraction of ALL the fibres of that fascicle):
      n_fibres, brain, muscle, frac_brain, frac_muscle
      with baseline: fired_prox (fibres the proximal pulse alone fires), blocked (of those, the ones
      that now deliver no AP to the muscle), block_ratio = blocked / fired_prox (the "anodic block"
      percentage of Tovbis et al. 2026), extra_brain (fibres that now send MORE APs to the brain than
      with the proximal pulse alone: extra sensations caused by the distal site), frac_extra_brain.
    """
    nb, nm = counts_at(ap, node_brain), counts_at(ap, node_muscle)
    if baseline is not None:
        b_nb, b_nm = counts_at(baseline, node_brain), counts_at(baseline, node_muscle)
        fired = (b_nb + b_nm) > 0
    rows = []
    for k in np.unique(fascicle):
        sel = fascicle == k
        n = int(sel.sum())
        row = {"fascicle": int(k), "n_fibres": n,
               "brain": int(np.sum(nb[sel] > 0)), "muscle": int(np.sum(nm[sel] > 0))}
        row["frac_brain"], row["frac_muscle"] = row["brain"] / n, row["muscle"] / n
        if baseline is not None:
            f = sel & fired
            row["fired_prox"] = int(f.sum())
            row["blocked"] = int(np.sum(f & (nm == 0)))
            row["block_ratio"] = row["blocked"] / row["fired_prox"] if row["fired_prox"] else np.nan
            row["extra_brain"] = int(np.sum(sel & (nb > b_nb)))
            row["frac_extra_brain"] = row["extra_brain"] / n
        rows.append(row)
    return pd.DataFrame(rows).set_index("fascicle")


# ----------------------------------------------------------------------------------------------------
# Space-time pictures of a fascicle
# ----------------------------------------------------------------------------------------------------
def volley_map(ap: list, fibres, n_nodes: int, frame_t_s: np.ndarray, window_s: float,
               n_norm: int | None = None) -> np.ndarray:
    """Fraction of the fibres carrying an AP at each node, at each frame time.

    A fibre "carries an AP at node n at time t" if its membrane crossed 0 mV at that node during the
    last window_s before t (an MRG AP stays above 0 mV for ~0.1-0.2 ms, so window_s ~ 0.1 ms draws
    each AP as a short bright segment that moves with the conduction velocity).

    ap         AP times of one scenario (all fibres)
    fibres     indices of the fibres to include (e.g. those of one fascicle)
    frame_t_s  (nT,) frame times [s]
    n_norm     denominator of the fraction (default: len(fibres)). Use the number of fibres of the
               fascicle so that the values are comparable with the 10 % criterion.
    Returns (nT, n_nodes) array in [0, 1].
    """
    fibres = list(fibres)
    n_norm = n_norm or max(1, len(fibres))
    out = np.zeros((len(frame_t_s), n_nodes))
    for i in fibres:
        for node, times in enumerate(ap[i]):
            if len(times) == 0:
                continue
            # frames in which at least one crossing of this node lies in (t - window, t]
            hit = np.zeros(len(frame_t_s), dtype=bool)
            for tc in times:
                hit |= (frame_t_s >= tc) & (frame_t_s < tc + window_s)
            out[hit, node] += 1
    return out / n_norm


def reach_profile(ap: list, fibres, n_nodes: int, t_from_s: float = -np.inf,
                  n_norm: int | None = None) -> np.ndarray:
    """For every node: fraction of the fibres with at least one AP at that node after t_from_s.

    Without block it is ~ the recruited fraction everywhere; an anodal block appears as a drop on
    the far side of the anode (the APs did not get past it).
    """
    fibres = list(fibres)
    n_norm = n_norm or max(1, len(fibres))
    hits = np.zeros(n_nodes)
    for i in fibres:
        for node, times in enumerate(ap[i]):
            if np.any(np.asarray(times) >= t_from_s):
                hits[node] += 1
    return hits / n_norm


def origin_nodes(ap: list) -> np.ndarray:
    """Node where each fibre's FIRST AP started (the earliest crossing of 0 mV), -1 if it never fires.

    Used as an artefact check (D-17): an AP that starts at one of the last nodes of the fibre was
    started by the cut end of the model fibre, not by the electrodes (real fibres do not end there).
    """
    out = np.full(len(ap), -1, dtype=int)
    for i, nodes in enumerate(ap):
        firsts = [(t[0], n) for n, t in enumerate(nodes) if len(t)]
        if firsts:
            out[i] = min(firsts)[1]
    return out


def started_at_end(ap: list, n_nodes: int, n_end: int = 3) -> np.ndarray:
    """True for the fibres whose first AP started within n_end nodes of either fibre end (D-17)."""
    o = origin_nodes(ap)
    return (o >= 0) & ((o < n_end) | (o > n_nodes - 1 - n_end))


def first_crossing_s(ap_fibre: list, node: int) -> float:
    """Time of the first AP at `node` for one fibre (NaN if none)."""
    t = ap_fibre[node]
    return float(t[0]) if len(t) else np.nan


# ----------------------------------------------------------------------------------------------------
# Safety (rough): Shannon criterion
# ----------------------------------------------------------------------------------------------------
def shannon_ok(charge_C: float, area_m2: float, k: float) -> tuple[bool, float, float]:
    """Shannon (1992) criterion log10(D) <= k - log10(Q), with Q in uC/phase, D in uC/cm^2/phase.

    Returns (within the limit?, Q [uC], D [uC/cm^2]). k = 1.75 is the value used for platinum by
    Tovbis et al. 2026 (citing Shannon 1992). The criterion was fitted on macro-electrodes; for small
    contacts like ours it is only indicative (docs/parameters.md).
    """
    Q_uC = charge_C * 1e6
    D_uC_cm2 = Q_uC / (area_m2 * 1e4)
    if Q_uC <= 0:
        return True, Q_uC, D_uC_cm2
    return bool(np.log10(D_uC_cm2) <= k - np.log10(Q_uC)), Q_uC, D_uC_cm2


# ----------------------------------------------------------------------------------------------------
# Export to the 3D viewer (playback of NEURON runs)
# ----------------------------------------------------------------------------------------------------
VIEWER_INDEX = "volleys.js"


def export_viewer_volley(out_dir, name: str, title: str, ap: list, *, fibre_x_m, fibre_y_m, fascicle,
                         z_nodes_m, node_brain: int, node_muscle: int, t_ref_s: float, channels: list,
                         target: int | None = None, z_sites_m=(), note: str = "") -> Path:
    """Write one simulated run for the viewer's Volley playback: <out_dir>/volley_<name>.js.

    The viewer (viewer/index.html, section "Volley") reads it with a <script> tag, like the field
    files, so the page still works when opened from disk. The file calls NV.addVolley({...}) with:
      fibres   x, y [mm] and fascicle of every fibre line (the fibres are straight lines along z)
      z_nodes  z [mm] of the nodes of Ranvier (the same for every fibre, nodeAlignment = proximalStation)
      events   the AP crossings as three parallel lists: fibre index, node index, time [ms] from t_ref
      channels the stimulator channels for the timeline: name, kind (cathode/anode), contacts
               {id: weight}, pulses [[t_on_ms, t_off_ms, amplitude_uA], ...] with times from t_ref
    and the index file <out_dir>/volleys.js (list of the runs present) is rewritten.

    out_dir   results/viewer/<patient>/ (the folder the viewer loads)
    name      file-name tag (letters, digits, _ . - only)
    ap        AP times of the run for ALL fibres (fibre_jobs.ap_times_job), times [s]
    t_ref_s   time origin of the playback [s] (the proximal onset)
    channels  [{"name", "kind", "contacts": {id: w}, "pulses": [(t_on_s, t_off_s, amp_A), ...]}, ...]
    """
    import json
    import re

    if not re.fullmatch(r"[A-Za-z0-9_.\-+]+", name):
        raise ValueError(f"volley name {name!r}: use letters, digits, _ . - + only")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # every AP crossing as one event (fibre, node, time in ms from t_ref_s): a compact flat list the viewer
    # indexes once (prepVolley in viewer/index.html)
    f_idx, n_idx, t_ms = [], [], []
    for i, nodes in enumerate(ap):
        for n, times in enumerate(nodes):
            for t in times:
                f_idx.append(i)
                n_idx.append(n)
                t_ms.append(round((float(t) - t_ref_s) * 1e3, 4))
    # geometry in mm, times in ms, currents in uA: the units of the viewer
    data = {
        "name": name, "title": title, "note": note,
        "target": None if target is None else int(target),
        "fibres": {"x": [round(float(v) * 1e3, 4) for v in fibre_x_m],
                   "y": [round(float(v) * 1e3, 4) for v in fibre_y_m],
                   "fascicle": [int(k) for k in fascicle]},
        "z_nodes": [round(float(z) * 1e3, 4) for z in z_nodes_m],
        "node_brain": int(node_brain), "node_muscle": int(node_muscle),
        "z_sites": [round(float(z) * 1e3, 4) for z in z_sites_m],
        "events": {"f": f_idx, "n": n_idx, "t": t_ms},
        "channels": [{"name": c["name"], "kind": c.get("kind", ""),
                      "contacts": {str(k): float(w) for k, w in c.get("contacts", {}).items()},
                      "pulses": [[round((on - t_ref_s) * 1e3, 5), round((off - t_ref_s) * 1e3, 5),
                                  round(a * 1e6, 4)] for on, off, a in c["pulses"]]} for c in channels],
    }
    path = out_dir / f"volley_{name}.js"
    path.write_text("NV.addVolley(" + json.dumps(data, separators=(",", ":")) + ");\n", encoding="utf-8")
    write_viewer_index(out_dir)
    return path


def write_viewer_index(out_dir) -> Path:
    """Rewrite <out_dir>/volleys.js with the names of the volley_<name>.js files present."""
    import json

    out_dir = Path(out_dir)
    names = sorted(p.stem[len("volley_"):] for p in out_dir.glob("volley_*.js"))
    path = out_dir / VIEWER_INDEX
    path.write_text("NV.setVolleyIndex(" + json.dumps(names) + ");\n", encoding="utf-8")
    return path


# ----------------------------------------------------------------------------------------------------
# Animation
# ----------------------------------------------------------------------------------------------------
def animate_profiles(path, curves: dict, frame_t_s: np.ndarray, z_nodes_m: np.ndarray, *,
                     ylabel: str, ylim: tuple, markers: list = (), title: str = "",
                     t_ref_s: float = 0.0, fps: int = 15, colours: dict | None = None,
                     styles: dict | None = None) -> Path:
    """Write an animated GIF of one or more profiles along the nerve.

    curves     {label: (nT, n_nodes) array}: one line per label, value at every node, per frame
    frame_t_s  (nT,) frame times [s]; the clock shows t - t_ref_s in ms
    markers    [(z_m, text), ...] vertical dashed lines (the two sites)
    Returns the path written. Needs Pillow (installed with matplotlib).
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    path = Path(path)
    zmm = z_nodes_m * 1e3
    fig, ax = plt.subplots(figsize=(9, 3.8))
    lines = {}
    for lab, arr in curves.items():
        kw = {}
        if colours and lab in colours:
            kw["color"] = colours[lab]
        if styles and lab in styles:
            kw["ls"] = styles[lab]
        (lines[lab],) = ax.plot(zmm, arr[0], lw=1.8, label=lab, **kw)
    for z, txt in markers:
        ax.axvline(z * 1e3, color="grey", ls="--", lw=0.8)
        ax.text(z * 1e3, ylim[1], txt, ha="center", va="bottom", fontsize=8, color="grey")
    clock = ax.text(0.01, 0.95, "", transform=ax.transAxes, va="top", fontsize=9)
    ax.set(xlim=(zmm[0], zmm[-1]), ylim=ylim, ylabel=ylabel,
           xlabel="z along the nerve [mm]   (brain <-  |  -> muscle)")
    ax.set_title(title, fontsize=9, pad=14)
    ax.legend(fontsize=7, loc="upper right", ncol=2)
    fig.tight_layout()

    def update(k):
        """Animation frame k: move every fascicle curve to its profile at frame_t_s[k], update the clock."""
        for lab, arr in curves.items():
            lines[lab].set_ydata(arr[k])
        clock.set_text(f"t = {(frame_t_s[k] - t_ref_s) * 1e3:+.2f} ms")
        return list(lines.values()) + [clock]

    anim = FuncAnimation(fig, update, frames=len(frame_t_s), blit=False)
    anim.save(str(path), writer=PillowWriter(fps=fps))
    plt.close(fig)
    return path
