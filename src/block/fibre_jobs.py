"""Batch runner for the MRG fibre simulations of the block study (serially or on several processes).

Why this module exists
    The notebooks of the block study (notebooks/block_anode_cathode.ipynb) need hundreds of fibre
    simulations. To use all the CPU cores they are sent to several worker processes. On Windows a
    worker process starts from scratch and can only run functions that live in a .py file it can
    import, NOT functions defined inside a notebook. So the few functions the workers need are here.

What a "job" is
    One simulation of one fibre, described by the tuple  (i_fibre, drives, t_stop_s):
      i_fibre   index of the fibre line in the profile file (0 ... nFibres-1)
      drives    list of (pattern, pulses), one entry per stimulator channel:
                  pattern = {contact id: weight}. The extracellular potential along the fibre for
                            1 A of this channel is  sum_c weight_c * V_c(z),  where V_c is the FEM
                            potential of contact c alone for 1 A (superposition: the FEM problem is
                            linear). The SIGN of the weight is the polarity of that contact:
                            cathodic (current drawn out of the tissue) < 0, anodic > 0.
                  pulses  = list of (t_on_s, t_off_s, amplitude_A) rectangular pulses, amplitude >= 0.
                            Pulses that overlap in time add up, so a ramp can be written as a
                            staircase of short pulses.
      t_stop_s  simulated time [s]
    The answer of a job (run_job) is (n_brain, n_muscle): the number of action potentials (APs) seen
    at the node detectInset_m from the proximal end of the fibre (z < 0, towards the brain) and at the
    node detectInset_m from the distal end (z > 0, towards the muscle). ap_times_job answers with the
    AP times at every node instead (population/volley analysis, notebooks/block_volley.ipynb).

One fibre object for all fibre lines
    Every process (the notebook's own and each worker) calls setup() once: it loads the FEM profiles
    and builds ONE MRG fibre in NEURON, which is then reused for every job. This is correct because
    all the fibre lines have the same node positions along z (a node exactly at the proximal station,
    b.mrg.nodeAlignment = "proximalStation" in config/block_params.m): the fibre lines differ only by
    the extracellular potential applied to them, which is set at the start of every job.

Units: SI everywhere (m, s, A, V), as in the rest of the project.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import numpy as np

from mrg_fibre import Channel, MRGFibre

# State of THIS process (filled by setup): profiles, the fibre object, the two detection nodes.
_S: dict = {}


def setup(h5_path, mrg: dict, diam_um: float, z_align_m: float, mech_dir) -> None:
    """Load the FEM profiles and build the MRG fibre in the current process.

    h5_path   results/block/<id>_fibre_profiles.h5 (written by scripts/export_block_profiles.m)
    mrg       the "mrg" block of the parameters (config/block_params.m -> <id>_block_params.json)
    diam_um   fibre diameter [um], one of the MRG_TABLE keys
    z_align_m z of one node of Ranvier [m] (the proximal station)
    mech_dir  folder with the compiled AXNODE mechanism (results/mrg/model)
    """
    with h5py.File(h5_path, "r") as f:
        _S["z_m"] = f["z_m"][:].ravel()                    # sample points along the nerve [m]
        _S["V"] = f["V_perA"][:]                           # (contacts, fibres, z) [V per A], float32
        _S["row"] = {int(c): k for k, c in enumerate(f["contact_id"][:].ravel())}
    _S["mrg"] = mrg
    fibre = MRGFibre(diam_um, mrg["fibreLength_m"], z_align_m, Path(mech_dir), mrg["temperature_C"])
    half_m = mrg["fibreLength_m"] / 2
    _S["fibre"] = fibre
    _S["node_brain"] = fibre.node_index(-half_m + mrg["detectInset_m"])   # proximal detection node
    _S["node_muscle"] = fibre.node_index(half_m - mrg["detectInset_m"])   # distal detection node


def pattern_profile(i_fibre: int, pattern: dict) -> np.ndarray:
    """Potential [V] along fibre i_fibre (at the z samples) for 1 A of the pattern current."""
    V, row = _S["V"], _S["row"]
    return sum(float(w) * V[row[int(c)], i_fibre, :].astype(float) for c, w in pattern.items())


def pattern_profiles_all(pattern: dict) -> np.ndarray:
    """Potential [V] along ALL fibres, shape (nFibres, nZ), for 1 A of the pattern current.

    Used for quick static screens (activating function on every fibre) without running NEURON.
    """
    V, row = _S["V"], _S["row"]
    return sum(float(w) * V[row[int(c)]].astype(float) for c, w in pattern.items())


def simulate(i_fibre: int, drives: list, t_stop_s: float, record: bool = False):
    """Run one fibre. Returns (n_brain, n_muscle, result dict of MRGFibre.run)."""
    m, fibre = _S["mrg"], _S["fibre"]
    channels = [Channel(pattern_profile(i_fibre, pat), list(pulses)) for pat, pulses in drives]
    res = fibre.run(channels, _S["z_m"], t_stop_s, m["dt_s"], m["vInit_V"], m["apThreshold_V"], record)
    n_brain = len(res["ap_times_s"][_S["node_brain"]])
    n_muscle = len(res["ap_times_s"][_S["node_muscle"]])
    return n_brain, n_muscle, res


def run_job(job) -> tuple[int, int]:
    """job = (i_fibre, drives, t_stop_s)  ->  (APs towards the brain, APs towards the muscle)."""
    i_fibre, drives, t_stop_s = job
    n_brain, n_muscle, _ = simulate(i_fibre, drives, t_stop_s)
    return n_brain, n_muscle


def ap_times_job(job) -> list:
    """job = (i_fibre, drives, t_stop_s)  ->  AP times at EVERY node of the fibre.

    Returns a list with one float32 array per node of Ranvier (node 0 = proximal end, towards the
    brain): the times [s] of the upward crossings of the AP threshold at that node. This is what the
    population ("volley") analysis of notebooks/block_volley.ipynb needs: from it one can tell where
    each AP travelled, where it stopped, and how many APs reached each end. It is small (a few numbers
    per node), so thousands of runs fit in memory. The counts of run_job are recovered as
    len(result[node_brain]) and len(result[node_muscle]).
    """
    i_fibre, drives, t_stop_s = job
    _, _, res = simulate(i_fibre, drives, t_stop_s)
    return [np.asarray(t, dtype=np.float32) for t in res["ap_times_s"]]


def profile_ap_times_job(job) -> list:
    """job = (channels, t_stop_s)  ->  AP times at every node, for potentials given EXPLICITLY.

    channels  list of (z_m, V_perA, pulses): the extracellular potential [V] along the fibre for +1 A,
              sampled at z_m [m] (any spacing, interpolated to the compartments), and its pulses as in
              a drive. Not tied to the FEM profiles: used for textbook checks of the fibre model (point
              electrodes in a homogeneous medium, notebooks/block_model_check.ipynb).
    The fibre is the shared one of this process (same diameter, length and node positions).
    """
    channels, t_stop_s = job
    m, fibre = _S["mrg"], _S["fibre"]
    z_ref = np.asarray(channels[0][0], dtype=float)
    chans = [Channel(np.interp(z_ref, np.asarray(z, float), np.asarray(v, float)), list(p))
             for z, v, p in channels]
    res = fibre.run(chans, z_ref, t_stop_s, m["dt_s"], m["vInit_V"], m["apThreshold_V"], False)
    return [np.asarray(t, dtype=np.float32) for t in res["ap_times_s"]]


def threshold_job(job) -> float:
    """Smallest amplitude [A] of one rectangular pulse on `pattern` that fires the fibre.

    job = (i_fibre, pattern, t_on_s, width_s, t_stop_s, start_A, max_A, rel_tol)
    "Fires" = at least one AP at either detection node. Search: start at start_A and double until the
    fibre fires (NaN if it still does not fire at max_A), then bisect down to rel_tol (relative).
    """
    i_fibre, pattern, t_on_s, width_s, t_stop_s, start_A, max_A, rel_tol = job

    def fires(amp_A):
        return sum(run_job((i_fibre, [(pattern, [(t_on_s, t_on_s + width_s, amp_A)])], t_stop_s))) > 0

    lo, hi = 0.0, start_A
    while not fires(hi):
        lo, hi = hi, 2 * hi
        if lo >= max_A:
            return float("nan")
    while hi - lo > rel_tol * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if fires(mid) else (mid, hi)
    return hi


class FibreRunner:
    """Runs lists of jobs, on n_workers processes (n_workers = 1: in the current process).

    The current process is always set up too, so that single runs with recordings (record) and the
    helper pattern_profile work in the notebook without going through the workers.
    """

    def __init__(self, h5_path, mrg: dict, diam_um: float, z_align_m: float, mech_dir, n_workers: int = 1):
        args = (str(h5_path), mrg, diam_um, z_align_m, str(mech_dir))
        setup(*args)
        self.n_workers = int(n_workers)
        self.pool = ProcessPoolExecutor(self.n_workers, initializer=setup, initargs=args) \
            if self.n_workers > 1 else None

    # the shared fibre object and the detection nodes, for plots and labels
    fibre = property(lambda self: _S["fibre"])
    node_brain = property(lambda self: _S["node_brain"])
    node_muscle = property(lambda self: _S["node_muscle"])

    def map(self, func, jobs: list, chunksize: int = 1) -> list:
        """[func(job) for job in jobs], in parallel if there are workers. func must live in this module."""
        if self.pool is None:
            return [func(j) for j in jobs]
        return list(self.pool.map(func, jobs, chunksize=chunksize))

    def run(self, jobs: list) -> list:
        """Run simulation jobs; returns [(n_brain, n_muscle), ...] in the order of the jobs."""
        return self.map(run_job, jobs, chunksize=4)

    def run_ap_times(self, jobs: list) -> list:
        """Run simulation jobs; returns the AP times at every node for each job (see ap_times_job)."""
        return self.map(ap_times_job, jobs, chunksize=4)

    def run_profiles(self, jobs: list) -> list:
        """Jobs with explicit potentials (see profile_ap_times_job); returns the AP times per node."""
        return self.map(profile_ap_times_job, jobs, chunksize=2)

    def record(self, i_fibre: int, drives: list, t_stop_s: float):
        """One run in the current process with the membrane voltage of every node recorded."""
        return simulate(i_fibre, drives, t_stop_s, record=True)

    def close(self) -> None:
        """Stop the worker processes (call it at the end of the notebook)."""
        if self.pool is not None:
            self.pool.shutdown()
            self.pool = None
