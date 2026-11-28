function validateQuicklookAnatomy(p) % Check editable nerve and fascicle inputs before plotting.
% Input p: quicklook parameter struct. If it carries a patient (p.geometry.patient), the patient is checked by validatePatient;
% otherwise the circular/elliptical layout p.geometry.fascicles_m (rows [x y a b angle_deg], m) is checked: positive sizes,
% every fascicle (with its perineurium, p.geometry.tPeri_m) inside the nerve ellipse, no two fascicles overlapping (tested
% conservatively with the larger semi-axis). No output: returns silently if valid, stops with an error otherwise.
if isfield(p.geometry,'patient') && ~isempty(p.geometry.patient), validatePatient(p.geometry.patient,p); return; end % A patient (virtual or measured) carries its own outlines and is validated as such.
fasc = p.geometry.fascicles_m; % Read the user-supplied [x y endoneurium radius] rows.
assert(size(fasc,2)==5 && size(fasc,1)>=1 && all(isfinite(fasc),'all'), ...
    'Fascicles must be rows [x y semi-axis-x semi-axis-y rotation-deg].'); % Reject malformed anatomy.
assert(all(fasc(:,3:4)>0,'all') && p.geometry.tPeri_m>0, 'Fascicle and perineurium sizes must be positive.'); % Require positive dimensions.
centres = fasc(:,1:2); % Extract the transverse fascicle centres.
distance = vecnorm(centres,2,2); % Measure each centre's distance from the nerve origin.
rays = centres ./ max(distance,eps); % Form the direction to each fascicle centre.
limit = nerveBoundaryRadius(p,rays); % Find the nerve contour in those directions.
assert(all(distance+max(fasc(:,3:4),[],2)+p.geometry.tPeri_m<limit), ...
    'A fascicle lies outside the supplied nerve contour.'); % Reject clearly infeasible placement.
for k = 1:size(fasc,1) % Check fascicle boundaries in the provisional circular representation.
    separation = vecnorm(centres(k+1:end,:)-centres(k,:),2,2); % Measure distances to later fascicles.
    minDistance = max(fasc(k,3:4))+max(fasc(k+1:end,3:4),[],2)+2*p.geometry.tPeri_m; % Include both perineurial sheaths.
    assert(all(separation>minDistance), 'The supplied fascicles overlap.'); % Require distinct regions.
end % Finish checking all fascicles.
end % Return after validating the anatomical input.
