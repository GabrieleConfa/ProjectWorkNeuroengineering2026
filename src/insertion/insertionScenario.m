function sc = insertionScenario(q) % One reproducible random scenario: irregular fascicles and one lead insertion per station.
%% Purpose
% This is the single source of truth for the production geometry in "insertion" mode: figures, mesh and viewer all
% read it. It draws random irregular fascicles (randomFascicleShapes) and, for every station, one random insertion
% of the PDMS lead from the entry wall to the opposite wall (simulateInsertion). Insertions that do not reach the far
% wall, or that leave a pad outside the nerve, are rejected and drawn again (fixed seeds, so the result is reproducible).
% The result is cached, because many functions ask for it.
persistent cache cacheKey % Memoization: the scenario is deterministic for the same inputs.
ip = insertionParams(); g = q.geometry; % Parameters and anatomy.
hasPatient = isfield(g,'patient') && ~isempty(g.patient); sig = []; if hasPatient, sig = cellfun(@(P)round(P(1:16:end,:),9),{g.patient.fascicles.outline},'UniformOutput',false); end
key = jsonencode({g.fascicles_m,g.nerveSemiAxes_m,g.nerveRotation_deg,g.tPeri_m,ip,q.electrode.reverseEntrySections,sig}); % Everything that changes the result.
if ~isempty(cache) && strcmp(cacheKey,key), sc = cache; return; end
%% 1. Random irregular fascicles (same centres and areas as the circular layout)
if hasPatient % Outlines of the patient (virtual or measured).
    pt = g.patient; polys0 = {pt.fascicles.outline}'; shapes = struct('polygons',{polys0},'nerve',pt.nerve.outline,'roundness',[pt.fascicles.roundness]','area_m2',[pt.fascicles.area_m2]');
else, shapes = randomFascicleShapes(q,ip,ip.seed); polys0 = shapes.polygons; end % Legacy: random shapes from the circular layout.
n = numel(polys0); % Undisturbed outlines.
raw = loadActiveSites(q.electrode.sitesFile,g.radius_m); sections = unique(raw.section)'; % Eight stations.
nerveIn = polyshape(shapes.nerve); % Nerve outline.
%% 2. One accepted insertion per station
stations = repmat(struct('section',0,'z_m',0,'u',[0 0],'polys',{{}},'shift',[],'path',[],'padXY',[],'padTan',[],'padTissue',[], ...
    'deflect',[],'displace',[],'pierce',[],'dirAngle_deg',0,'attempts',0),1,numel(sections)); % Preallocate.
for iS = 1:numel(sections)
    sec = sections(iS); rows = find(raw.section==sec & raw.type=="intraneural"); sites = raw(rows,:); % The four internal pads of this station.
    [~,i] = max(sites.r_m); u = sites.centre_m(i,1:2)/sites.r_m(i); % Nominal insertion axis (as in the original lead planner).
    if ismember(sec,q.electrode.reverseEntrySections), u = -u; end % Same entry side as the approved lead trajectories.
    padS = sites.centre_m(:,1:2)*u'; tipS = nerveBoundaryRadius(q,u); entryS = -nerveBoundaryRadius(q,-u); % Pad positions, far wall and entry wall.
    st = struct('u',u,'padS',padS,'tipS',tipS,'entryS',entryS,'nerve',shapes.nerve,'section',sec,'padHalf',sites.sideXY_m(1)/2); % Station description.
    accepted = false;
    for attempt = 1:3000 % Draw until a valid insertion is found.
        rs = RandStream('mt19937ar','Seed',ip.seed+1000*iS+attempt); % Reproducible stream per attempt.
        o = simulateInsertion(q,ip,polys0,st,rs); % One random insertion.
        if ~o.completed || any(isnan(o.padXY(:,1))), continue; end % The lead must reach the far wall and carry all four pads.
        okPads = true; % Pads must lie inside the nerve and not overlap each other.
        for k = 1:size(o.padXY,1)
            e = o.padXY(k,:)+[-1;1]*st.padHalf*o.padTan(k,:); % The two ends of the active face.
            okPads = okPads && all(isinterior(nerveIn,e(:,1),e(:,2)));
        end
        if size(o.padXY,1)>1, dmin = min(pdist1(o.padXY)); okPads = okPads && dmin > 2*st.padHalf*1.05; end % Pads must not overlap.
        if ~okPads, continue; end
        accepted = true; break
    end
    assert(accepted,'No valid lead insertion at station %+d after 3000 random draws: the section is probably too crowded near the nerve wall for a lead that touches both walls. Try another patient or relax the constraints.',sec);
    stations(iS) = struct('section',sec,'z_m',sites.z_m(1),'u',u,'polys',{o.polys},'shift',o.shift,'path',anchorPath(o.path,ip.stepLength_m,0), ...
        'padXY',o.padXY,'padTan',o.padTan,'padTissue',o.padTissue,'deflect',o.deflect,'displace',o.displace,'pierce',o.pierce, ...
        'dirAngle_deg',atan2d(o.direction*[-u(2);u(1)],o.direction*u'),'attempts',attempt); % Store the accepted insertion.
end
%% 3. Output
sc = struct('ip',ip,'polys0',{polys0},'nerve',shapes.nerve,'roundness',shapes.roundness,'area_m2',shapes.area_m2,'stations',stations, ...
    'bump_m',0.4e-3,'leadWidth_m',ip.tipDiameter_m,'leadThickness_m',0.10e-3,'sleeveThickness_m',0.15e-3); % Axial extent of a displacement, PDMS lead and sleeve sizes.
cache = sc; cacheKey = key; % Remember it.
end
function d = pdist1(P) % Distances between all pairs of points (no toolbox needed).
n = size(P,1); d = []; for i = 1:n-1, for j = i+1:n, d(end+1) = norm(P(i,:)-P(j,:)); end, end %#ok<AGROW>
end
function P = anchorPath(path,ds,extra) % Decimate the path (keeping the exact anchored end) and prolong the entry end so the lead crosses the entry wall.
seg = vecnorm(diff(path),2,2); a = [0;cumsum(seg)]; keep = [true;diff(floor(a/(10*ds)))>0]; keep(end) = true; P = path(keep,:); % About 100 um spacing.
t0 = P(2,:)-P(1,:); t0 = t0/norm(t0); % Entry direction.
if extra>0, P = [P(1,:)-extra*t0; P]; end % Optional prolongation of the entry end (not used: both ends stop against the wall).
end
