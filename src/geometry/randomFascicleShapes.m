function shapes = randomFascicleShapes(q, ip, seed) % Replace the circular fascicles by random irregular smooth shapes.
%% Purpose
% Fascicles are neither circles nor ellipses (tutor). Each shape keeps the area and the centre of the circular
% layout in q.geometry.fascicles_m, but gets a random elongation (roundness drawn from the histology statistics) and a
% random smooth irregularity (radial Fourier harmonics). Shapes never overlap and stay inside the nerve.
%% Inputs / outputs
%   q       quicklook parameter struct: q.geometry.fascicles_m (rows [x y a b angle_deg], m), tPeri_m [m], nerve ellipse
%   ip      insertionParams(): roundnessMean, roundnessSD (Pusnik 2025), harmonicAmplitude, fascicleGap_m [m]
%   seed    seed of the random stream (same seed -> same shapes)
%   shapes  struct: polygons (n x 1 cell of 128 x 2 outlines [m]), nerve (720 x 2 outline [m]), roundness (n x 1,
%           4A/(pi*major^2)), area_m2 (n x 1, equal to the circular areas), attempts (draws used per fascicle)
%% 1. Inputs and containers
g = q.geometry; F = g.fascicles_m; n = size(F,1); % Circular layout: [x y a b angle] per fascicle.
nPts = 128; th = (0:nPts-1)'*2*pi/nPts; % Vertices per outline (dense enough for distances and meshing).
rs = RandStream('mt19937ar','Seed',seed); % Reproducible random numbers.
area = pi*F(:,3).*F(:,4); [~,order] = sort(area,'descend'); % Place the largest fascicle first.
gap = 2*g.tPeri_m + ip.fascicleGap_m; % Minimum distance between outlines (two perineuria plus a small gap).
nerve = nerveOutline(q,720); % Nerve outline polygon.
inner = polybuffer(polyshape(nerve),-(g.tPeri_m+0.035e-3)); % Allowed region for fascicle outlines.
polys = cell(n,1); roundness = zeros(n,1); attempts = zeros(n,1); % Outputs.
%% 2. Random shape for each fascicle, rejected until it fits
for k = order(:)' % Visit fascicles from the largest.
    A = area(k); placed = find(~cellfun(@isempty,polys)); % Target area and shapes already placed.
    for t = 1:2500 % Rejection sampling with progressive relaxation.
        if t<=500, mu = ip.roundnessMean; sd = ip.roundnessSD; amp = ip.harmonicAmplitude; % Statistics from the literature.
        elseif t<=1500, mu = 0.8; sd = 0.05; amp = 0.5*ip.harmonicAmplitude; % Fall back to rounder shapes.
        else, mu = 0.95; sd = 0.02; amp = 0.01; end % Last resort: almost a circle.
        rd = min(0.95,max(0.5,mu+sd*randn(rs))); % Target roundness (4A/(pi*major^2) = 1/axis ratio for an ellipse).
        AR = 1/rd; a = sqrt(A*AR/pi); b = sqrt(A/(pi*AR)); % Semi-axes of the base ellipse with exactly the target area.
        harm = zeros(nPts,1); % Radial irregularity.
        for m = 3:5, harm = harm + amp/(m-2)*(2*rand(rs)-1)*cos(m*th+2*pi*rand(rs)); end % Decaying random harmonics, m = 3..5.
        P = [a*cos(th), b*sin(th)].*(1+harm); % Base ellipse deformed radially.
        P = P*sqrt(A/polyarea(P(:,1),P(:,2))); % Restore the exact area.
        phi = 2*pi*rand(rs); P = P*[cos(phi) sin(phi);-sin(phi) cos(phi)]; % Random orientation.
        jitter = ((t>100)*0.04e-3+(t>600)*0.06e-3+(t>1200)*0.10e-3)*randn(rs,1,2); % Centre jitter, larger after many failed attempts (dense layouts).
        P = P+F(k,1:2)+jitter; % Move to the layout centre.
        feret = max(max(sqrt((P(:,1)-P(:,1)').^2+(P(:,2)-P(:,2)').^2))); % Maximum Feret diameter.
        r = 4*A/(pi*feret^2); % Roundness as defined in the histology paper.
        if t<=500 && abs(r-ip.roundnessMean)>2*ip.roundnessSD, continue; end % Keep the literature statistics (within 2 SD) while the first attempts last.
        if ~all(inpolygon(P(:,1),P(:,2),inner.Vertices(:,1),inner.Vertices(:,2))), continue; end % Must stay inside the nerve.
        ok = true; % Check the distance to the fascicles already placed.
        for j = placed(:)'
            if polyDist(P,polys{j})<gap || any(inpolygon(P(:,1),P(:,2),polys{j}(:,1),polys{j}(:,2))) || any(inpolygon(polys{j}(:,1),polys{j}(:,2),P(:,1),P(:,2))), ok = false; break; end
        end
        if ok, polys{k} = P; roundness(k) = r; attempts(k) = t; break; end % Accept this shape.
    end
    if isempty(polys{k}) % Very dense layout: use the exact circle of the layout (always valid), so the patient can still be built.
        polys{k} = F(k,1:2)+F(k,3)*[cos(th) sin(th)]; roundness(k) = 1; attempts(k) = 0;
    end
end
%% 3. Output
shapes = struct('polygons',{polys},'roundness',roundness,'area_m2',area,'attempts',attempts,'nerve',nerve,'seed',seed); % Shared by the insertion study and the figures.
end
function d = polyDist(A,B) % Minimum distance between two closed polylines (vertex-to-segment, both directions).
d = min(ptSeg(A,B),ptSeg(B,A)); % Symmetric.
end
function d = ptSeg(P,Q) % Minimum distance from the points P to the closed polyline Q.
Q2 = circshift(Q,-1,1); S = Q2-Q; L2 = sum(S.^2,2)'; % Segment vectors and squared lengths.
px = P(:,1); py = P(:,2); % Points as columns.
t = ((px-Q(:,1)').*S(:,1)'+(py-Q(:,2)').*S(:,2)')./L2; t = min(1,max(0,t)); % Projection parameter on each segment.
dx = px-(Q(:,1)'+t.*S(:,1)'); dy = py-(Q(:,2)'+t.*S(:,2)'); % Offsets to the closest point of each segment.
d = min(sqrt(dx.^2+dy.^2),[],'all'); % Smallest distance.
end
function xy = nerveOutline(q,n) % Nerve outline of the preview anatomy as a polygon.
t = (0:n-1)'*2*pi/n; ray = [cos(t) sin(t)]; r = nerveBoundaryRadius(q,ray); xy = r(:).*ray; % Polar radius of the oval along each ray.
end
