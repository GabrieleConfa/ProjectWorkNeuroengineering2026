function [fig, model] = run_3d_quicklook(mode, p) % Show the 3D nerve and its end-on device projection.
arguments % Set defaults while allowing the anatomical inputs to be replaced later.
    mode (1,1) string {mustBeMember(mode, ["deflect","penetrate"])} = "deflect" % Select the insertion scenario.
    p (1,1) struct = quicklook3DParams() % Load the editable geometry defaults.
end % Finish describing the inputs.

%% 1. Load contacts and assemble the geometry in metres
validateQuicklookAnatomy(p); % Check anatomical dimensions and placement before any geometry work.
sites = loadActiveSites(p.electrode.sitesFile, p.geometry.radius_m); % Load the tutor's finite contact dimensions and orientations.
cuffSites = sites(sites.type == "cuff", :); % Select the two external contacts per station.
stripSites = sites(sites.type == "intraneural", :); % Select the four internal pads per station.
model.nerve = buildQuicklookNerve(p, cuffSites); % Flatten the nerve against the rigid external contacts.
ids = unique(stripSites.section, 'sorted'); % Enumerate the eight axial stations.
for k = 1:numel(ids) % Build a transverse strip independently at each station.
    strip = stripSites(stripSites.section == ids(k), :); % Extract the internal pads for this station.
    cuff = cuffSites(cuffSites.section == ids(k), :); % Extract the corresponding opposite cuff pair.
    model.strips(k) = planQuicklookStrip(p, strip, cuff, mode); % Place the strip and record its geometric overlaps.
    model.puncture(k) = assessFasciclePuncture(p, strip); % Screen the straight encounter with the provisional bending inputs.
end % Finish assembling the device.
model.parameters = p; % Preserve the input parameters with the output.
model.mode = mode; % Preserve the interaction assumption for later comparison.

%% 2. Draw matching 3D and end-on views
fig = plotQuicklookNerve(model); % Return a MATLAB figure that can be rotated and zoomed.
end % Return the figure and geometry for inspection or later mesh work.
