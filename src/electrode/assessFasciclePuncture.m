function assessment = assessFasciclePuncture(p, sites) % Screen whether the free strip can deliver a reported puncture load.
%% 1. Read dimensions and compute transverse bending resistance
width = p.electrode.stripWidth_m; % Read the strip width from one editable parameter field.
thickness = p.electrode.stripThickness_m; % Read the strip thickness from one editable parameter field.
modulus = p.electrode.stripYoungModulus_Pa; % Read the assumed strip elastic modulus.
length = p.electrode.unsupportedLength_m; % Read the unsupported length that can buckle.
threshold = p.electrode.perineuriumPunctureForce_N; % Read the experimental comparison load.
assert(all([width,thickness,modulus,length,threshold] > 0), 'Mechanical inputs must be positive.'); % Reject unspecified or invalid inputs.
I_m4 = width*thickness^3/12; % Compute the weak-axis second moment for a rectangular strip.
Pcritical_N = pi^2*modulus*I_m4/(4*length^2); % Compute ideal fixed-free Euler buckling load.
%% 2. Find which fascicles the original straight insertion would meet
[~,first] = max(sites.r_m); % Select an outer contact to set the straight insertion direction.
u = sites.centre_m(first,1:2)/sites.r_m(first); % Recover the tutor's nominal transverse insertion axis.
fasc = p.geometry.fascicles_m; % Read fascicle centres and endoneurium radii.
along_m = fasc(:,1:2)*u'; % Locate each fascicle along the insertion axis.
normal_m = fasc(:,1:2) - along_m*u; % Measure the transverse offset from that axis.
outer_m = max(fasc(:,3:4),[],2) + p.geometry.tPeri_m; % Include the perineurium around each fascicle.
intersects = vecnorm(normal_m,2,2) <= outer_m + width/2; % Include finite strip width in the straight-path collision.
entry_s = -nerveBoundaryRadius(p,-u); % Keep encounters on the physical one-sided insertion only.
overhang = 0.05e-3; % Read legacy cases with the provisional default.
if isfield(p.electrode,'tipOverhang_m'), overhang=p.electrode.tipOverhang_m; end % Accept the measured tip extension when supplied.
tip_s = max(sites.centre_m(:,1:2)*u'+sites.sideXY_m/2)+overhang; % Stop beyond the fourth active face.
halfChord = sqrt(max(0,(outer_m+width/2).^2-sum(normal_m.^2,2))); % Find both limits of contact with the conservative envelope.
intersects = intersects & along_m+halfChord>=entry_s & along_m-halfChord<=tip_s; % Exclude obstacles beyond the device.
ids = find(intersects); % Report the fascicles met by the undeviated strip.
normalAngle_deg = nan(size(ids)); % Allocate incidence angles for the encountered fascicles.
for k = 1:numel(ids) % Compute a normal incidence estimate at each first contact.
    j = ids(k); % Read this fascicle's index.
    radial = sqrt(max(0,(outer_m(j)+width/2)^2 - norm(normal_m(j))^2)); % Find the axial entry offset at the expanded boundary.
    entry = (along_m(j)-radial)*u; % Locate the first contact point along the insertion line.
    n = (entry-fasc(j,1:2)); % Compute the boundary normal at that entry point.
    normalAngle_deg(k) = acosd(min(1,abs(dot(u,n)/norm(n)))); % Measure incidence from the surface normal.
end % Finish measuring the encounter angles.
%% 3. Compare mechanical capacity with the measured probe load
assessment.section = sites.section(1); % Preserve the electrode station identifier.
assessment.fascicleIds = ids; % Report straight-path encounters for interpretation.
assessment.normalAngle_deg = normalAngle_deg; % Report angle without claiming an angle-dependent puncture law.
assessment.Pcritical_N = Pcritical_N; % Return the ideal unsupported-strip buckling load.
assessment.comparisonForce_N = threshold; % Return the experimental proxy load.
assessment.loadRatio = Pcritical_N/threshold; % Compare device buckling with the probe's measured puncture force.
assessment.inference = "uncertain"; % Start with an unresolved outcome because no calibrated contact law exists.
if isempty(ids) % Handle a path that never meets a fascicle.
    assessment.inference = "no straight-path encounter"; % Report a purely geometric absence of contact.
elseif Pcritical_N < threshold % Check whether an unsupported strip would buckle before the proxy puncture load.
    assessment.inference = "deflection expected for unsupported strip"; % Report a conditional mechanics screen.
end % Keep larger loads unresolved until the insertion force and tip geometry are supplied.
end % Return all intermediate quantities for review.
