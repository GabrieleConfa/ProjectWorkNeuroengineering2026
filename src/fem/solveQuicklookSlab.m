function result = solveQuicklookSlab(p, meshPreview, cuffSites) % Solve a local two-pad electrical check on the contact slab.
%% 1. Assign the three tissue conductivities from the source material table
X = meshPreview.mesh.Nodes; % Read tetrahedral coordinates in metres.
T = meshPreview.mesh.Elements(1:4,:); % Use linear tetrahedra for the shared P1 assembler.
bt = meshPreview.regionTable; % Read each cell's anatomical shape membership.
nE = size(T,2); % Count tetrahedra for vectorized conductivity assignment.
sigma = zeros(3,nE); % Allocate transverse and longitudinal conductivity per element.
roles = zeros(size(bt,1),1); % Keep a transparent tissue label for each geometry cell.
for cellId = 1:size(bt,1) % Assign every source geometry cell independently.
    members = find(bt(cellId,:)); % Read nerve, perineurium and endoneurium membership.
    if any(mod(members,2)==1 & members>1) % Find the odd-index endoneurium shapes.
        role = 1; % Select the endoneurial conductivity.
    elseif any(mod(members,2)==0 & members>1) % Find the even-index perineurial shapes.
        role = 2; % Select the perineurial conductivity.
    else % Treat tissue between fascicles as epineurium.
        role = 3; % Select the epineurial conductivity.
    end % Finish anatomical classification of this cell.
    elem = findElements(meshPreview.mesh,"region","Cell",cellId); % Locate all tetrahedra in the cell.
    sigma(:,elem) = repmat([p.materials.sigmaT_Spm(role);p.materials.sigmaT_Spm(role);p.materials.sigmaL_Spm(role)],1,numel(elem)); % Preserve longitudinal anisotropy.
    roles(cellId) = role; % Retain the cell tissue type for inspection.
end % Finish material assignment.
assert(all(sigma(:)>0), 'Every tetrahedron needs a positive tissue conductivity.'); % Reject incomplete cell mapping.
%% 2. Set the two finite cuff faces as ideal surface electrodes
assert(height(cuffSites)==2, 'The local solve needs the opposing cuff pad pair.'); % Keep a well-defined two-terminal check.
active = padNodes(X,p,cuffSites(1,:)); % Find nerve nodes on the first finite pad footprint.
returnNodes = padNodes(X,p,cuffSites(2,:)); % Find nerve nodes on the opposing pad footprint.
assert(~isempty(active) && ~isempty(returnNodes) && isempty(intersect(active,returnNodes)), ...
    'The finite pad footprints did not map to distinct mesh boundary nodes.'); % Reject an invalid electrode mapping.
%% 3. Assemble, solve once and scale the linear field to one microampere
K = assembleConductionP1(X,T,sigma); % Reuse the colleague's validated tetrahedral conduction assembly.
u = zeros(size(X,2),1); % Set the return pad to zero and prepare a one-volt drive.
u(active) = 1; % Hold the first pad equipotential without modelling its bulk material.
fixed = unique([active(:);returnNodes(:)]); % Identify both Dirichlet electrode boundaries.
free = setdiff((1:size(X,2))',fixed); % Leave all other tissue surfaces naturally insulating.
u(free) = K(free,free)\(-K(free,fixed)*u(fixed)); % Solve the stationary volume-conductor equation.
drive_A = sum(K(active,:)*u); % Calculate current implied by the one-volt solution.
assert(isfinite(drive_A) && drive_A>0, 'The two-terminal current must be positive.'); % Reject a disconnected solve.
V = u*(1e-6/drive_A); % Scale the potential to a one-microampere input.
return_A = sum(K(returnNodes,:)*V); % Check the integrated current at the return pad.
assert(abs(return_A+1e-6)<1e-8, 'Current does not balance at the return pad.'); % Enforce current conservation.
%% 4. Return the field and its explicit modelling limits
result.nodes_m = X; % Return the original spatial mesh nodes.
result.tets = T; % Return the original tetrahedra.
result.V_per_uA_V = V; % Return potential for one microampere, in volts.
result.driveCurrentAtOneVolt_A = drive_A; % Expose the model conductance for checking.
result.returnCurrent_A = return_A; % Expose the signed current balance.
result.activeNodes = active; % Preserve the finite active electrode footprint.
result.returnNodes = returnNodes; % Preserve the finite return electrode footprint.
result.cellRole = roles; % Preserve endoneurium/perineurium/epineurium labels.
result.scope = "local slab only"; % Prevent interpretation as the full-length stimulation model.
end % Finish the local electrical verification.

function nodes = padNodes(X,p,site) % Identify boundary mesh nodes under one rigid contact face.
n = [cosd(site.theta_deg),sind(site.theta_deg)]; % Read the outward radial direction of this cuff contact.
t = [-n(2),n(1)]; % Read the transverse direction across the pad.
plane = cuffContactPlane(p,site); % Compute the same contact plane used to build the mesh.
normal_m = n*X(1:2,:); % Project every mesh node onto the pad normal.
lateral_m = t*X(1:2,:); % Measure node position across the contact width.
tol_m = max(1e-6,p.meshPreview.hMin_m/20); % Allow the polygon approximation of a curved surface.
atPlane = abs(normal_m-plane)<tol_m; % Keep only nodes lying on the flat contact surface.
inWidth = abs(lateral_m)<=site.sideXY_m/2+tol_m; % Enforce the measured finite transverse width.
inLength = abs(X(3,:)-site.z_m)<=site.sideZ_m/2+tol_m; % Enforce the measured finite axial length.
nodes = find(atPlane & inWidth & inLength); % Return all active surface nodes on this pad.
end % Finish mapping one electrode face to mesh nodes.
