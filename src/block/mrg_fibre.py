"""MRG myelinated fibre driven by the FEM extracellular potential (conduction-block study).

Purpose
    The static FEM gives the potential V(z) along a fibre; whether an action potential (AP) starts,
    travels, collides with another AP or is blocked under an anode is a question of membrane
    dynamics. This module builds the MRG double-cable fibre model (McIntyre, Richardson & Grill 2002,
    J Neurophysiol 87:995, ModelDB 3810) in NEURON and applies to it any combination of contact
    currents with their own time courses.

    The morphology and the electrical parameters are those of MRGaxon.hoc (ModelDB 3810), rebuilt
    here in Python so that the fibre length (number of nodes) and the diameter are parameters; only
    the nodal membrane mechanism (AXNODE.mod) is compiled, in results/mrg/model/ (see docs/block.md).

Units
    Inputs and outputs in SI (m, s, A, V), as in the rest of the project. NEURON works in um, ms, mV:
    the conversion is done here and nowhere else.

Main objects
    MRGFibre(diam_um, length_m, z_align_m, mech_dir)   one fibre along z, centred on z = 0
    Channel(profile_V_perA, pulses)                    one stimulator output: the potential along the
                                                       fibre for +1 A (any weighted sum of contacts)
                                                       and its rectangular current pulses
    MRGFibre.run(channels, t_stop_s, ...)              simulate, return membrane voltage at the nodes
                                                       and the AP crossing times
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from neuron import h

# Geometry table of MRGaxon.hoc (proc dependent_var), one row per fibre diameter fiberD [um]:
# g ratio, axon diameter, node diameter, MYSA diameter, FLUT diameter [um], internode deltax [um],
# FLUT length paralength2 [um], number of myelin lamellae nl.
MRG_TABLE = {
    5.7:  (0.605, 3.4, 1.9, 1.9, 3.4, 500, 35, 80),
    7.3:  (0.630, 4.6, 2.4, 2.4, 4.6, 750, 38, 100),
    8.7:  (0.661, 5.8, 2.8, 2.8, 5.8, 1000, 40, 110),
    10.0: (0.690, 6.9, 3.3, 3.3, 6.9, 1150, 46, 120),
    11.5: (0.700, 8.1, 3.7, 3.7, 8.1, 1250, 50, 130),
    12.8: (0.719, 9.2, 4.2, 4.2, 9.2, 1350, 54, 135),
    14.0: (0.739, 10.4, 4.7, 4.7, 10.4, 1400, 56, 140),
    15.0: (0.767, 11.5, 5.0, 5.0, 11.5, 1450, 58, 145),
    16.0: (0.791, 12.7, 5.5, 5.5, 12.7, 1500, 60, 150),
}
# Constants of MRGaxon.hoc (proc model_globels)
NODE_LENGTH_UM = 1.0
MYSA_LENGTH_UM = 3.0          # paralength1
SPACE_P1_UM, SPACE_P2_UM, SPACE_I_UM = 0.002, 0.004, 0.004   # periaxonal space widths
RHOA_OHM_UM = 0.7e6           # axoplasmic resistivity
MYELIN_CM = 0.1               # uF/cm2 per lamella membrane
MYELIN_GM = 0.001             # S/cm2 per lamella membrane

_LOADED: set[str] = set()


def load_mechanism(mech_dir: Path) -> None:
    """Load the compiled AXNODE mechanism once (nrnivmodl output in mech_dir)."""
    mech_dir = Path(mech_dir).resolve()
    if str(mech_dir) in _LOADED:
        return
    if sys.platform.startswith("win"):
        dll = mech_dir / "nrnmech.dll"
        if not dll.exists():
            raise FileNotFoundError(f"{dll} not found: compile AXNODE.mod with nrnivmodl (docs/block.md)")
        h.nrn_load_dll(str(dll))
    else:
        import neuron
        neuron.load_mechanisms(str(mech_dir))
    h.load_file("stdrun.hoc")
    _LOADED.add(str(mech_dir))


@dataclass
class Channel:
    """One stimulator output.

    profile_V_perA : potential [V] at the fibre sample points z_m for +1 A (weighted sum of contacts;
                     the weights carry the polarity: cathodic = negative weight)
    pulses         : list of (t_on_s, t_off_s, amplitude_A) rectangular pulses
    """
    profile_V_perA: np.ndarray
    pulses: list[tuple[float, float, float]] = field(default_factory=list)


class MRGFibre:
    """One MRG fibre along z, nodes at z_align_m + k * internode, spanning |z| <= length_m / 2."""

    def __init__(self, diam_um: float, length_m: float, z_align_m: float, mech_dir: Path,
                 temperature_C: float = 37.0):
        load_mechanism(mech_dir)
        if diam_um not in MRG_TABLE:
            raise ValueError(f"MRG has no parameter set for D = {diam_um} um; allowed {sorted(MRG_TABLE)}")
        g, axonD, nodeD, paraD1, paraD2, deltax, paralength2, nl = MRG_TABLE[diam_um]
        self.diam_um, self.internode_m = diam_um, deltax * 1e-6
        h.celsius = temperature_C

        # nodes inside |z| <= L/2, one of them exactly at z_align_m
        half_um, za_um = length_m * 1e6 / 2, z_align_m * 1e6
        k_back = int(np.floor((za_um + half_um) / deltax))
        z_first = za_um - k_back * deltax
        n_nodes = int(np.floor((half_um - z_first) / deltax)) + 1
        if n_nodes < 5:
            raise ValueError("fibre too short")
        self.n_nodes = n_nodes

        # resistances of the periaxonal space (MOhm/cm) and internodal segment length, as in the hoc
        pi = np.pi
        Rpn0 = (RHOA_OHM_UM * .01) / (pi * (((nodeD / 2) + SPACE_P1_UM) ** 2 - (nodeD / 2) ** 2))
        Rpn1 = (RHOA_OHM_UM * .01) / (pi * (((paraD1 / 2) + SPACE_P1_UM) ** 2 - (paraD1 / 2) ** 2))
        Rpn2 = (RHOA_OHM_UM * .01) / (pi * (((paraD2 / 2) + SPACE_P2_UM) ** 2 - (paraD2 / 2) ** 2))
        Rpx = (RHOA_OHM_UM * .01) / (pi * (((axonD / 2) + SPACE_I_UM) ** 2 - (axonD / 2) ** 2))
        interlength = (deltax - NODE_LENGTH_UM - 2 * MYSA_LENGTH_UM - 2 * paralength2) / 6
        xg_my, xc_my = MYELIN_GM / (nl * 2), MYELIN_CM / (nl * 2)
        D = diam_um

        def make(kind, i):
            s = h.Section(name=f"{kind}{i}")
            s.nseg = 1
            if kind == "node":
                s.diam, s.L, s.Ra, s.cm = nodeD, NODE_LENGTH_UM, RHOA_OHM_UM / 10000, 2
                s.insert("axnode")
                s.insert("extracellular")
                s(0.5).xraxial[0], s(0.5).xg[0], s(0.5).xc[0] = Rpn0, 1e10, 0
                return s
            if kind == "MYSA":
                dr, L, gp, rp = paraD1, MYSA_LENGTH_UM, 0.001 * paraD1 / D, Rpn1
            elif kind == "FLUT":
                dr, L, gp, rp = paraD2, paralength2, 0.0001 * paraD2 / D, Rpn2
            else:  # STIN
                dr, L, gp, rp = axonD, interlength, 0.0001 * axonD / D, Rpx
            s.diam, s.L = D, L
            s.Ra = RHOA_OHM_UM * (1 / (dr / D) ** 2) / 10000
            s.cm = 2 * dr / D
            s.insert("pas")
            s(0.5).g_pas, s(0.5).e_pas = gp, -80
            s.insert("extracellular")
            s(0.5).xraxial[0], s(0.5).xg[0], s(0.5).xc[0] = rp, xg_my, xc_my
            return s

        # chain: node, MYSA, FLUT, 6 STIN, FLUT, MYSA, node, ... (connections as in MRGaxon.hoc)
        self.nodes = [make("node", i) for i in range(n_nodes)]
        chain, k = [self.nodes[0]], 0
        for i in range(n_nodes - 1):
            seq = [make("MYSA", 2 * i), make("FLUT", 2 * i)] + [make("STIN", 6 * i + j) for j in range(6)] \
                + [make("FLUT", 2 * i + 1), make("MYSA", 2 * i + 1)]
            prev = self.nodes[i]
            for s in seq:
                s.connect(prev(1), 0)
                prev = s
            self.nodes[i + 1].connect(prev(1), 0)
            chain += seq + [self.nodes[i + 1]]
        self.sections = chain
        lengths = np.array([s.L for s in chain])
        starts = np.concatenate([[0.0], np.cumsum(lengths)[:-1]])
        centre_um = starts + lengths / 2 - NODE_LENGTH_UM / 2 + z_first     # node 0 centred at z_first
        self.z_sec_m = centre_um * 1e-6
        self.z_node_m = np.array([self.z_sec_m[chain.index(n)] for n in self.nodes])
        self._segs = [s(0.5) for s in chain]

    # ------------------------------------------------------------------------------------------------
    def node_index(self, z_m: float) -> int:
        """Index of the node closest to z_m."""
        return int(np.argmin(np.abs(self.z_node_m - z_m)))

    def run(self, channels: list[Channel], z_profile_m: np.ndarray, t_stop_s: float, dt_s: float = 5e-6,
            v_init_V: float = -80e-3, ap_threshold_V: float = 0.0, record_all: bool = False) -> dict:
        """Simulate the fibre under the given channels.

        The extracellular potential at each compartment is sum_c I_c(t) * profile_c(z), with the
        profiles interpolated linearly from the samples z_profile_m to the compartment centres.
        Returns a dict with
          t_s           (nT,) time [s] (only if record_all)
          v_V           (nNodes, nT) membrane voltage at the nodes [V] (only if record_all)
          ap_times_s    list (nNodes) of arrays: times of the upward crossings of ap_threshold_V
        """
        prof = np.array([np.interp(self.z_sec_m, z_profile_m, c.profile_V_perA) for c in channels])
        # breakpoints of the piecewise-constant currents
        edges = {0.0, t_stop_s}
        for c in channels:
            for t_on, t_off, _ in c.pulses:
                edges.update([t_on, t_off])
        edges = sorted(t for t in edges if 0 <= t <= t_stop_s)

        h.dt = dt_s * 1e3
        tvec = h.Vector().record(h._ref_t)
        vvecs = [h.Vector().record(n(0.5)._ref_v) for n in self.nodes]
        for seg in self._segs:
            seg.e_extracellular = 0.0
        h.finitialize(v_init_V * 1e3)
        for t0, t1 in zip(edges[:-1], edges[1:]):
            tm = 0.5 * (t0 + t1)
            I = np.array([sum(a for on, off, a in c.pulses if on <= tm < off) for c in channels])
            ve_mV = (I @ prof) * 1e3 if len(channels) else np.zeros(len(self._segs))
            for seg, v in zip(self._segs, ve_mV):
                seg.e_extracellular = v
            h.continuerun(t1 * 1e3)
        for seg in self._segs:
            seg.e_extracellular = 0.0

        t = tvec.as_numpy().copy() * 1e-3
        V = np.array([v.as_numpy() for v in vvecs]) * 1e-3
        up = (V[:, :-1] < ap_threshold_V) & (V[:, 1:] >= ap_threshold_V)
        ap_times = [t[1:][up[i]] for i in range(V.shape[0])]
        res = {"ap_times_s": ap_times}
        if record_all:
            res["t_s"], res["v_V"] = t, V
        return res


def count_aps(res: dict, node: int, t_from_s: float = 0.0) -> int:
    """Number of APs (upward crossings) recorded at a node after t_from_s."""
    return int(np.sum(res["ap_times_s"][node] >= t_from_s))
