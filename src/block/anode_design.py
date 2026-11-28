"""Design of the distal blocking configuration by optimisation on the FEM solution (D-19).

Why this module exists
    The mirror anode (same contact as the proximal cathode, D-16/D-18) never blocked the volley of the
    target fascicle: a 0.2 mm contact hyperpolarises ONE node of Ranvier (the AP jumps over it) and, as
    its current rises, some other spot depolarises first and fires the fibres (virtual cathodes, sleeve
    ends, return contacts; docs/block.md section 6). This module searches, among ALL the contacts, the
    configuration (which contacts, signed currents) that best does what a block needs:

      * hyperpolarise K consecutive nodes of every fibre that the proximal pulse fires in the target
        fascicle, as strongly as possible ("block depth" h);
      * keep the depolarisation of EVERY node of EVERY fibre (all fascicles) below one common limit d
        (otherwise the distal site fires fibres itself: extra sensation, twitch);
      * do not hyperpolarise the target fibres on the brain side of the stimulation site (the sensation
        must pass), here at most a fraction of h;
      * use at most N contacts (the "minimum electrode": a mixed-integer constraint).

    The quantity it maximises is the ratio h / d. Because the extracellular potential is linear in the
    contact currents (superposition, D-07), the activating function (AF) of any configuration is
    A @ I, with A the AF of each contact for 1 A; the problem is then a mixed-integer linear program
    (scipy.optimize.milp, HiGHS). The ratio is scale free: it does not need the (unknown) absolute AF
    thresholds of excitation and block. Which ratio is enough, and at which current, is decided by the
    MRG fibre model in NEURON (notebooks/block_anode_design.ipynb). The AF is only a proxy (it ignores
    the pulse shape and how wide a depolarised region is), so the optimiser proposes and NEURON checks.

Units: SI (A, m, V/m^2) unless the name says otherwise.
"""
from __future__ import annotations

import h5py
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from mrg_fibre import MRG_TABLE


def node_positions(diam_um: float, length_m: float, z_align_m: float) -> np.ndarray:
    """z [m] of the nodes of Ranvier of the MRG fibre (same rule as mrg_fibre.MRGFibre)."""
    dx_um = MRG_TABLE[diam_um][5]
    half_um, za_um = length_m * 1e6 / 2, z_align_m * 1e6
    k_back = int(np.floor((za_um + half_um) / dx_um))
    z_first = za_um - k_back * dx_um
    n_nodes = int(np.floor((half_um - z_first) / dx_um)) + 1
    return (z_first + dx_um * np.arange(n_nodes)) * 1e-6


def unit_af(h5_path, z_nodes_m: np.ndarray) -> dict:
    """Activating function at the nodes for 1 A in each contact alone.

    Returns a dict with
      A         (nContacts, nFibres, nNodes) [V/m^2 per A]: second difference of the FEM potential over
                one internode at the interior nodes; at the two end nodes the sealed-end form
                (V_{1} - V_{0}) / dx^2 (the fibre ends are part of the model, D-17)
      contact_id, fascicle, x_m, y_m  (from the profile file)
    """
    with h5py.File(h5_path, "r") as f:
        z = f["z_m"][:].ravel()
        V = f["V_perA"][:].astype(float)               # (C, F, Z) [V per A]
        out = {"contact_id": f["contact_id"][:].ravel().astype(int),
               "fascicle": f["fascicle"][:].ravel().astype(int),
               "x_m": f["x_m"][:].ravel(), "y_m": f["y_m"][:].ravel()}
    dx = float(np.mean(np.diff(z_nodes_m)))
    Ve = np.empty(V.shape[:2] + (len(z_nodes_m),))
    for c in range(V.shape[0]):
        for i in range(V.shape[1]):
            Ve[c, i] = np.interp(z_nodes_m, z, V[c, i])
    A = np.empty_like(Ve)
    A[:, :, 1:-1] = (Ve[:, :, :-2] - 2 * Ve[:, :, 1:-1] + Ve[:, :, 2:]) / dx ** 2
    A[:, :, 0] = (Ve[:, :, 1] - Ve[:, :, 0]) / dx ** 2
    A[:, :, -1] = (Ve[:, :, -2] - Ve[:, :, -1]) / dx ** 2
    out["A"] = A
    return out


def design(A: np.ndarray, contact_id: np.ndarray, *, allowed, target_fibres, block_nodes,
           protect_nodes=(), protect_frac: float = 0.3, n_max: int | None = None,
           balanced: bool = False, time_limit_s: float = 60.0) -> dict | None:
    """Best configuration for one choice of block nodes (and at most n_max contacts).

    A              unit AF (nContacts, nFibres, nNodes) [V/m^2 per A]
    contact_id     id of each row of A
    allowed        ids of the contacts the optimiser may use
    target_fibres  fibre indices whose APs must be blocked (the fibres the proximal pulse fires)
    block_nodes    node indices to hyperpolarise (K consecutive nodes)
    protect_nodes  node indices on the brain side of the stimulation site: on the target fibres the
                   hyperpolarisation there must stay below protect_frac x h (the sensation must pass)
    n_max          at most this many active contacts (None: no limit, a plain LP)
    balanced       True: sum of the currents = 0 (no net current to the far ground)

    Normalisation: the depolarisation of every node of every fibre is <= 1 (arbitrary unit), so the
    optimum h is directly the ratio block depth / worst depolarisation. Returns a dict
      ratio     h (hyperpolarisation of the block nodes on the worst target fibre / worst depolarisation)
      weights   {contact id: weight}, scaled so that the largest |weight| is 1 (sign = polarity:
                anode > 0, cathode < 0); the stimulator amplitude multiplies these weights
      af_dep_per_A  worst depolarisation [V/m^2] for 1 A of the largest contact (to convert ratios
                    into currents: the NEURON sweep starts from here)
      n_active  number of contacts used
    or None if the solver failed.
    """
    allowed = [c for c in allowed if c in set(contact_id.tolist())]
    rows_c = np.array([list(contact_id).index(c) for c in allowed])
    Ac = A[rows_c]                                     # (nc, F, N)
    nc = len(rows_c)
    scale = np.abs(Ac).max()                           # work with O(1) numbers
    Ac = Ac / scale
    use_int = n_max is not None
    nv = nc + 1 + (nc if use_int else 0)               # [I (nc), h, z (nc)]
    D_all = Ac.reshape(nc, -1).T                       # (every fibre x node, nc): depolarisation rows
    # (1) "depolarisation <= 1" has one row per node of every fibre (~25 000), but only a few are ever
    # binding. Constraint generation: start with the rows each contact pushes hardest, solve, add the rows
    # the solution violates, solve again (exact at convergence, much faster than all rows at once).
    work = set()
    for j in range(nc):
        work.update(np.argsort(-D_all[:, j])[:200].tolist())
    # per-contact current bound: 5 x the current that contact could carry alone without exceeding the
    # depolarisation limit (in either polarity). A tight bound makes the mixed-integer search fast; the
    # factor 5 leaves room for contacts that cancel each other's depolarisation.
    pos = 1.0 / np.maximum(D_all.max(0), 1e-9)
    neg = 1.0 / np.maximum((-D_all).max(0), 1e-9)
    i_bound = 5.0 * np.maximum(pos, neg)
    # phase 1 on the continuous relaxation (fast), phase 2 with the integer variables (if n_max is given):
    # the rows collected in phase 1 are usually enough, so phase 2 needs only one or two integer solves
    phases = [(False, nc + 1, 25)] + ([(True, nv, 6)] if use_int else [])
    for integer, nvar, max_iter in phases:
        for _ in range(max_iter):
            res, rows_used = _solve(Ac, D_all, sorted(work), nc, nvar, integer, target_fibres, block_nodes,
                                    protect_nodes, protect_frac, balanced, n_max, time_limit_s, i_bound)
            if res.x is None:
                return None
            # constraint generation: add the depolarisation rows the solution violates (up to 3000 per
            # round, worst first) and solve again, until no row of the full set is violated
            af = D_all @ res.x[:nc]
            viol = np.flatnonzero(af > 1.0 + 1e-6)
            new = [r for r in viol[np.argsort(-af[viol])] if r not in work][:3000]
            if not new:
                break
            work.update(new)
    # result: contact currents normalised to the largest one (= weight 1), contacts below 1e-4 dropped
    I = res.x[:nc]
    h = float(res.x[nc])
    big = np.abs(I).max()
    if big <= 0:
        return None
    keep = np.abs(I) > 1e-4 * big
    weights = {int(allowed[j]): float(I[j] / big) for j in range(nc) if keep[j]}
    wvec = np.array([weights.get(c, 0.0) for c in allowed])
    dep_per_A = float(np.max(np.tensordot(wvec, A[rows_c], axes=1)))
    return {"ratio": h, "weights": weights, "af_dep_per_A": dep_per_A, "n_active": int(keep.sum()),
            "block_nodes": [int(m) for m in block_nodes], "optimal": bool(res.success),
            "rows_used": rows_used}


def _solve(Ac, D_all, work, nc, nv, use_int, target_fibres, block_nodes, protect_nodes, protect_frac,
           balanced, n_max, time_limit_s, i_bound):
    """One MILP/LP solve with the depolarisation rows in `work` (see design)."""
    rows, lo, hi = [], [], []

    def add(M, l, u):
        """Append a block of linear constraints  l <= M @ x <= u  (rows of M, bounds per row)."""
        rows.append(M); lo.append(l); hi.append(u)

    D = D_all[work]
    M = np.zeros((D.shape[0], nv)); M[:, :nc] = D
    add(M, np.full(D.shape[0], -np.inf), np.ones(D.shape[0]))
    # (2) block nodes of the target fibres: AF + h <= 0
    H = np.array([Ac[:, i, m] for i in target_fibres for m in block_nodes])
    M = np.zeros((len(H), nv)); M[:, :nc] = H; M[:, nc] = 1.0
    add(M, np.full(len(H), -np.inf), np.zeros(len(H)))
    # (3) brain side of the target fibres: AF + protect_frac h >= 0
    if len(protect_nodes):
        B = np.array([Ac[:, i, m] for i in target_fibres for m in protect_nodes])
        M = np.zeros((len(B), nv)); M[:, :nc] = B; M[:, nc] = protect_frac
        add(M, np.zeros(len(B)), np.full(len(B), np.inf))
    # (4) balanced: sum I = 0
    if balanced:
        M = np.zeros((1, nv)); M[0, :nc] = 1.0
        add(M, [0.0], [0.0])
    lb = np.r_[-i_bound, 0.0]
    ub = np.r_[i_bound, np.inf]
    integ = np.zeros(nv)
    if use_int:                                        # |I_c| <= i_max z_c, sum z <= n_max
        E = np.zeros((2 * nc, nv))
        for j in range(nc):
            E[2 * j, j] = 1.0;  E[2 * j, nc + 1 + j] = -i_bound[j]
            E[2 * j + 1, j] = -1.0; E[2 * j + 1, nc + 1 + j] = -i_bound[j]
        add(E, np.full(2 * nc, -np.inf), np.zeros(2 * nc))
        S = np.zeros((1, nv)); S[0, nc + 1:] = 1.0
        add(S, [-np.inf], [float(n_max)])
        integ[nc + 1:] = 1
        lb = np.r_[lb, np.zeros(nc)]; ub = np.r_[ub, np.ones(nc)]
    c_obj = np.zeros(nv); c_obj[nc] = -1.0             # maximise h
    res = milp(c_obj, constraints=LinearConstraint(np.vstack(rows), np.concatenate(lo), np.concatenate(hi)),
               integrality=integ, bounds=Bounds(lb, ub), options={"time_limit": time_limit_s, "disp": False})
    return res, len(work)


def af_of(A: np.ndarray, contact_id: np.ndarray, weights: dict) -> np.ndarray:
    """AF (nFibres, nNodes) [V/m^2] of a configuration, for 1 A times its weights."""
    idx = {int(c): k for k, c in enumerate(contact_id)}
    return sum(w * A[idx[int(c)]] for c, w in weights.items())
