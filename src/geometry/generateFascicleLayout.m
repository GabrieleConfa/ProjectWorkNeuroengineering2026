function [R_m, fasc, info] = generateFascicleLayout(gen, peri, footprints_m)
%GENERATEFASCICLELAYOUT Random, literature-based cross-section: circular nerve + circular fascicles.
%   [R_m, fasc, info] = generateFascicleLayout(gen, peri, footprints_m)
%
%   Draws the fascicles of one virtual patient. The caller sets the random stream (rng(seed))
%   beforehand, so the layout is reproducible from the seed alone.
%
%   Biology built in (sources in config/pw6_3d_params.m and docs/parameters.md):
%     - nerve area, fascicle count, total fascicle area and the one large fascicle: Pusnik 2025;
%     - the other fascicles have log-normal equivalent diameters with the coefficient of variation
%       of the upper-arm median nerve (Brill & Tyler 2017), scaled to the total area;
%     - positions spread homogeneously over the section, regardless of size (Delgado-Martinez 2016):
%       uniform random start, then overlapping fascicles are pushed apart ("relaxation");
%     - number of inner fascicles (not on the convex hull of the centres) ~ 0.58 (N - 6.08)
%       (Brill & Tyler 2017).
%   Fascicles keep a band of epineurium between them and from the nerve surface (gen.minGap*).
%
%   Shank contacts (gen.contactRule):
%     "pierce"    the shank goes through the tissue, as a transversal electrode does: a contact may
%                 lie in a fascicle, in the epineurium, or across a fascicle boundary. Only nearly
%                 tangent configurations (a fascicle boundary closer than gen.contactClearance_m to
%                 a contact edge or corner without crossing it cleanly) are rejected: they would
%                 create sliver elements in the mesh.
%     "pushAside" the tissue is pushed aside by the rigid contact: every contact lies entirely in one
%                 fascicle or entirely in the epineurium. With the fascicles extruded along z, all 16
%                 shank footprints apply to the same section, which forces an electrode-shaped layout
%                 (big fascicle on the axis, the others centred on contacts): kept only for comparison.
%
%   Inputs
%     gen           p.geometry.generator (see config/pw6_3d_params.m)
%     peri          p.perineurium (realThickness_m = anatomical thickness, used for the gaps)
%     footprints_m  cell array, one 4x2 matrix [x y] [m] per shank-contact footprint in the x-y plane
%
%   Outputs
%     R_m   nerve radius [m]
%     fasc  table, one row per fascicle, largest first: x_m, y_m, rEndo_m (endoneurium radius = where
%           the perineurium thin layer is), tPeri_m (anatomical perineurium thickness)
%     info  struct: restarts, iterations, nInner, targetInner, failures

arguments
    gen (1,1) struct
    peri (1,1) struct
    footprints_m cell
end

R_m = sqrt(gen.nerveArea_m2 / pi);
N   = gen.nFascicles;
tP  = peri.realThickness_m;
clr = gen.contactClearance_m;
targetInner = round(0.58 * (N - 6.08));
sigmaLn = sqrt(log(1 + gen.diameterCV^2));            % log-normal shape for the given CV
nFp = numel(footprints_m);
fpCentre = cell2mat(cellfun(@(c) mean(c, 1), footprints_m(:), 'UniformOutput', false));
fpHalfDiag = cellfun(@(c) max(vecnorm(c - mean(c, 1), 2, 2)), footprints_m(:));
switch gen.contactRule
    case "pierce",    fpOK = @(c, r, fp) pierceOK(c, r, fp, clr);
    case "pushAside", fpOK = @(c, r, fp) pushAsideOK(c, r, r + tP, fp, clr);
    otherwise, error('generateFascicleLayout:rule', 'Unknown contactRule "%s".', gen.contactRule);
end

nFail = zeros(1, 4);                                  % sizes, relaxation, final check, inner rule
for restart = 1:gen.maxRestarts
    % --- sizes ---------------------------------------------------------------------------------
    d = exp(sigmaLn * randn(N - 1, 1));               % relative diameters (median 1)
    a = d.^2 / sum(d.^2) * (gen.totalFascicleArea_m2 - gen.largestFascicleArea_m2);
    rEndo = sqrt([gen.largestFascicleArea_m2; sort(a, 'descend')] / pi);
    if any(rEndo(2:end) >= rEndo(1)) || any(rEndo < gen.minEndoRadius_m)
        nFail(1) = nFail(1) + 1;   continue           % redraw the sizes
    end
    rOut = rEndo + tP;                                % anatomical outer radius, for the gaps

    % --- positions: uniform start, then relaxation ---------------------------------------------
    rMax = R_m - gen.minGapNerve_m - rOut;
    if any(rMax <= 0), nFail(1) = nFail(1) + 1; continue; end
    rho = rMax .* sqrt(rand(N, 1));   th = 2*pi*rand(N, 1);
    c = [rho .* cos(th), rho .* sin(th)];
    for it = 1:gen.maxRelaxIterations
        moved = false;
        for i = 1:N                                   % pairs: push apart, the larger moves less
            for j = i+1:N
                v = c(i,:) - c(j,:);   dist = norm(v);
                need = rOut(i) + rOut(j) + gen.minGapFascicles_m;
                if dist < need
                    if dist == 0, v = [cos(it) sin(it)]; dist = 1; end
                    u = v / dist;   wI = rOut(j)^2 / (rOut(i)^2 + rOut(j)^2);
                    step = 1.02 * (need - dist);
                    c(i,:) = c(i,:) + u * step * wI;
                    c(j,:) = c(j,:) - u * step * (1 - wI);
                    moved = true;
                end
            end
        end
        for i = 1:N                                   % stay inside the nerve
            lim = R_m - gen.minGapNerve_m - rOut(i);
            if norm(c(i,:)) > lim, c(i,:) = c(i,:) * lim / norm(c(i,:)); moved = true; end
        end
        for i = 1:N                                   % contacts
            shift = [0 0];
            for k = 1:nFp
                if fpOK(c(i,:), rEndo(i), footprints_m{k}), continue; end
                v = c(i,:) - fpCentre(k,:);   dist = max(norm(v), eps);   u = v / dist;
                if gen.contactRule == "pierce"        % small step away from the bad configuration
                    shift = shift + u * 2 * clr * sign(rEndo(i) - dist + eps);
                else
                    inLim  = rEndo(i) - clr - fpHalfDiag(k);
                    outLim = rOut(i) + clr + fpHalfDiag(k);
                    if inLim > 0 && dist < rEndo(i), shift = shift + u * (0.97 * inLim - dist);
                    else,                            shift = shift + u * (1.03 * outLim - dist);
                    end
                end
                moved = true;
            end
            c(i,:) = c(i,:) + shift;
        end
        if moved && mod(it, 300) == 0                 % stuck: small random kick
            c = c + 0.05 * mean(rOut) * randn(N, 2);
        end
        if ~moved, break; end
    end
    if moved, nFail(2) = nFail(2) + 1; continue; end
    if ~layoutOK(c, rEndo, rOut, R_m, gen, footprints_m, fpOK), nFail(3) = nFail(3) + 1; continue; end

    % --- biology check: inner fascicles --------------------------------------------------------
    hull = unique(convhull(c(:,1), c(:,2)));
    nInner = N - numel(hull);
    if gen.useInnerFascicleRule && abs(nInner - targetInner) > gen.innerFascicleTolerance
        nFail(4) = nFail(4) + 1;   continue
    end
    fasc = table(c(:,1), c(:,2), rEndo, repmat(tP, N, 1), ...
                 'VariableNames', {'x_m', 'y_m', 'rEndo_m', 'tPeri_m'});
    info = struct('restarts', restart, 'iterations', it, 'nInner', nInner, ...
                  'targetInner', targetInner, 'failures', nFail);
    return
end
error('generateFascicleLayout:failed', ...
      ['No valid layout after %d restarts (failures: sizes %d, relaxation %d, final check %d, ' ...
       'inner-fascicle rule %d): relax the clearances.'], gen.maxRestarts, nFail(1), nFail(2), nFail(3), nFail(4));
end

function ok = pierceOK(c, r, corners, clr)
%PIERCEOK The circle (centre c, radius r) is not nearly tangent to the rectangle: its distance from
% every corner and from the closest point of every edge differs from r by at least clr.
dCorner = vecnorm(corners - c, 2, 2);
dEdge = zeros(4, 1);
for k = 1:4
    a = corners(k,:);   b = corners(mod(k, 4) + 1, :);
    t = max(0, min(1, dot(c - a, b - a) / dot(b - a, b - a)));
    dEdge(k) = norm(c - (a + t * (b - a)));
end
ok = all(abs([dCorner; dEdge] - r) >= clr);
end

function ok = pushAsideOK(c, rIn, rOut, corners, clr)
%PUSHASIDEOK The rectangle is inside the circle rIn - clr or outside the circle rOut + clr.
dMax = max(vecnorm(corners - c, 2, 2));
ok = dMax <= rIn - clr || pointRectDistance(c, corners) >= rOut + clr;
end

function d = pointRectDistance(p, corners)
%POINTRECTDISTANCE Distance from point p (1x2) to the convex polygon "corners" (4x2); 0 if inside.
if inpolygon(p(1), p(2), corners(:,1), corners(:,2)), d = 0; return; end
d = inf;
for k = 1:4
    a = corners(k,:);   b = corners(mod(k, 4) + 1, :);
    t = max(0, min(1, dot(p - a, b - a) / dot(b - a, b - a)));
    d = min(d, norm(p - (a + t * (b - a))));
end
end

function ok = layoutOK(c, rEndo, rOut, R, gen, footprints, fpOK)
%LAYOUTOK Exact final check of all the constraints.
N = size(c, 1);   ok = true;
for i = 1:N
    ok = ok && norm(c(i,:)) + rOut(i) <= R - gen.minGapNerve_m + 1e-12;
    for j = i+1:N
        ok = ok && norm(c(i,:) - c(j,:)) >= rOut(i) + rOut(j) + gen.minGapFascicles_m - 1e-12;
    end
    for k = 1:numel(footprints)
        ok = ok && fpOK(c(i,:), rEndo(i), footprints{k});
    end
end
end
