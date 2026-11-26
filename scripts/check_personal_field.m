function checks = check_personal_field(matFile) % Verify the source labels, contact locations and conservation after the run.
%% 1. Check correspondence between solved fields and geometry
S = load(matFile,'out'); % Read the numerical result saved before graphical export.
o = S.out; % Access the full solution object.
assert(o.patient.id=="Cuff60" && isfield(o.params,'personalGeometry'),'Expected the personal 3D solution.'); % Reject reference data.
q = o.params.personalGeometry; % Read the anatomical input actually saved by the solver.
assert(abs((max(o.mesh.nodes_m(3,:))-min(o.mesh.nodes_m(3,:)))-q.domain.length_m)<1e-9,'Unexpected domain length.'); % Verify both ends of the real mesh.
assert(isequal(sort(o.solvedIds(:)),sort(o.sites.id)),'Missing contact solution.'); % Require the requested set of 48 columns.
assert(all(isfinite(o.V_per_A),'all') && all(o.results.R_kOhm>0),'Invalid contact potential.'); % Verify finite values and positive input resistance.
raw = loadActiveSites(q.electrode.sitesFile,q.geometry.radius_m); % Recover the nominal angles and straight pad positions.
largestError_m = 0; % Accumulate the largest preview-to-FEM centre displacement.
for station = unique(raw.section)' % Audit each axial level separately.
    internal = raw.section==station & raw.type=="intraneural"; % Locate the four penetrating pads.
    external = raw.section==station & raw.type=="cuff"; % Locate the two rigid pads.
    strip = planQuicklookStrip(q,raw(internal,:),raw(external,:),"deflect"); % Recreate the intended device from editable settings.
    if ~isfield(q.geometry,'scenario') || isempty(q.geometry.scenario), assert(isempty(strip.envelopeFascicles),'Intersection detected in the deflection scenario.'); end % Legacy planner: require clearance. In insertion mode a lead may cross a fascicle (outcome of the insertion).
    xy = squeeze(mean(strip.padFaces_m,1))'; % Extract the four active face centres.
    largestError_m = max(largestError_m,max(vecnorm(xy-o.sites.centre_m(internal,:),2,2))); % Compare to actual FEM input coordinates.
end % Finish all eight sections.
assert(largestError_m<1e-10,'Preview and solved contacts disagree.'); % Reject a stale or mismatched geometry.
%% 2. Report the numerical errors on the tested mesh
checks = table(height(o.results),size(o.mesh.tets,2),largestError_m,max(o.results.conservErr_pct),max(o.results.spread_pct),'VariableNames',{'solvedContacts','tetrahedra','centreError_m','maxCurrentError_pct','maxContactSpread_pct'}); % Summarize reproducible diagnostics.
assert(checks.maxCurrentError_pct<0.01,'Current balance error exceeds the numerical validation tolerance.'); % Require input and return agreement.
if isfield(o,'cuffAudit'), plotPersonalMesh(o,fileparts(matFile)); exportPersonalMeshView(o,fileparts(matFile)); end % Export the actual mesh only after verification.
disp(checks); % Show actual measured errors for the simulation report.
end % Return the checked metrics to the user or caller.
