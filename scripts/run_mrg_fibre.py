"""Couple a MATLAB axial FEM profile to the published ModelDB 3810 motor fibre."""
import argparse  # Accept the exported field and adjustable pulse parameters.
from pathlib import Path  # Resolve the local ModelDB checkout.
import numpy as np  # Interpolate the field and record numerical output.
from neuron import h  # Use the compiled original MRG membrane mechanism.

# 1. Load the published 10 µm motor fibre (37 °C, 21 nodes).
parser = argparse.ArgumentParser(description=__doc__)  # Describe this single-fibre experiment.
parser.add_argument("profile", type=Path)  # Two columns from export_axon_profile.
parser.add_argument("--current-uA", type=float, default=100)  # Exploratory phase amplitude.
parser.add_argument("--phase-ms", type=float, default=0.05)  # Exploratory phase duration.
args = parser.parse_args()  # Read the requested stimulation scenario.
assert args.current_uA > 0 and args.phase_ms > 0  # Require positive phase magnitudes.
folder = Path(__file__).resolve().parents[1] / "results/mrg"  # Locate the local model and compiled mechanism.
h.load_file("stdrun.hoc")  # Load fixed-step recording and simulation procedures.
h.load_file(str(folder / "model/MRGaxon.hoc"))  # Instantiate the published double-cable morphology.
h.stim.amp = 0  # Disable the original intracellular current injection.
h.distance(0, 0.5, sec=h.node[0])  # Measure arc length from the first node centre.
centre = h.distance(0.5, sec=h.node[10])  # Centre the fibre at the active contact's axial station.
sections = list(h.allsec())  # Include nodes, paranodes and internodes.
z = np.array([h.distance(0.5, sec=s)-centre for s in sections])  # Positions in micrometres.
profile = np.loadtxt(args.profile, delimiter=",")  # Read axial position and extracellular volts per ampere.
assert np.isfinite(profile).all() and np.all(np.diff(profile[:,0])>0)  # Verify the full input profile.
assert z.min() >= profile[0,0] and z.max() <= profile[-1,0]  # Reject extrapolation outside the solved nerve.
transfer = np.interp(z, profile[:,0], profile[:,1])  # Sample the FEM at every electrical compartment.

# 2. Apply a balanced cathodic-first pulse and a zero-current control.
h.dt = min(0.005,args.phase_ms/20)  # Resolve each phase with at least twenty integration steps.
indices = (10,13,16,19)  # Record the central node and three progressively more distal nodes.
def simulate(amplitude):
    times = h.Vector().record(h._ref_t)  # Record elapsed milliseconds.
    traces = [h.Vector().record(h.node[i](0.5)._ref_v) for i in range(int(h.axonnodes))]  # Record membrane voltage at each monitored node.
    h.finitialize(-80)  # Start from the published resting voltage.
    for phase, end in [(0,1),(-1,1+args.phase_ms),(1,1+2*args.phase_ms),(0,5+2*args.phase_ms)]:
        for sec, gain in zip(sections,transfer):  # Apply the local potential at each compartment.
            sec(0.5).e_extracellular = gain*amplitude*phase*1e-3  # Convert volts per ampere times microamperes to mV.
        h.continuerun(end)  # Integrate until the next phase transition.
    return np.column_stack([times.as_numpy().copy()]+[v.as_numpy().copy() for v in traces])  # Preserve this run.

control = simulate(0)  # Check for spontaneous spikes without stimulation.
assert np.max(control[:,1:])<0  # Require the unstimulated control to stay below the spike detection level.
trace = simulate(args.current_uA)  # Apply the requested extracellular field waveform.
base = args.profile.with_name(f"{args.profile.stem}_{args.current_uA:g}uA_{args.phase_ms:g}ms")  # Identify the patient, contact and fascicle in the result.
np.savetxt(str(base)+"_traces.csv",trace,delimiter=",",header="t_ms,"+",".join(f"node{i}_mV" for i in range(int(h.axonnodes))),comments="")  # Save every membrane trace.
first_times = np.full(int(h.axonnodes), np.inf)  # Keep every nodal crossing time to detect end initiation.
for node in range(int(h.axonnodes)):
    col = node+1  # Column zero holds simulation time.  # Report independent times of upward zero crossings.
    rises = np.flatnonzero((trace[:-1,col]<0)&(trace[1:,col]>=0))  # Detect candidate action potentials.
    if rises.size: first_times[node] = trace[rises[0]+1,0]  # Save the earliest upward crossing.
    if node in indices: print(node,trace[rises+1,0].tolist())  # Distal crossings at increasing times indicate propagation.

if np.isfinite(first_times).any():  # Locate the earliest spike before interpreting distal latency.
    first = np.flatnonzero(first_times==first_times.min())  # Include ties at the simulation time resolution.
    print("First crossing at node(s)",first.tolist(),"of",len(first_times),"from 0")  # Flag boundary initiation explicitly.
    if np.any((first<=1)|(first>=len(first_times)-2)): print("End initiation: threshold is sensitive to fibre truncation.")  # Avoid interpreting an endpoint spike as local recruitment.

# 3. Save the propagation graph for comparison with the static FEM maps.
import matplotlib  # Select an offscreen backend for batch execution.
matplotlib.use("Agg")  # Save the figure without opening a second GUI.
import matplotlib.pyplot as plt  # Plot voltages with units.
fig, ax = plt.subplots()  # Create one axis for all four recording locations.
for node in indices:  # Keep the same node order as the CSV file.
    ax.plot(trace[:,0],trace[:,node+1],label=f"Node {node}")  # Plot the recorded transmembrane voltage.
ax.set(xlabel="Time (ms)",ylabel="Membrane potential (mV)",title=f"{base.name} | MRG 10 µm | {args.current_uA:g} µA, {args.phase_ms:g} ms/phase")  # Label the exploratory parameters.
ax.legend()  # Identify the recording sites.
fig.tight_layout()  # Keep physical units visible.
fig.savefig(str(base)+"_propagation.png",dpi=160)  # Export the evidence for this single fibre.
