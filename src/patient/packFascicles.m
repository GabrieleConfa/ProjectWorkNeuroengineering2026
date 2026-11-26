function F = packFascicles(q, rs) % Random layout of circular fascicles inside the nerve by random start and repulsion (relaxation).
%% Purpose
% Random sequential placement jams at about 55 percent coverage, which is where the aggregate anatomy of the mid-arm nerve lies
% (fascicles cover 54 percent of the section). Instead, fascicles start at random positions and are pushed apart, and away from
% the nerve surface, until nothing overlaps. Areas: the largest fascicle and the total area are fixed by the histology
% (Pusnik 2025), the other nine share the rest with random unequal weights. Output rows: [x y r r 0] in metres.
g = q.geometry; spec = g.pw6; n = spec.nFascicles; tP = g.tPeri_m;
w = sort(0.5+rand(rs,n-1,1),'descend'); areas = [spec.areaLargest_m2;(spec.areaTotal_m2-spec.areaLargest_m2)*w/sum(w)]; r = sqrt(areas/pi); % Areas of the fascicles.
th = (0:359)'*pi/180; ray = [cos(th) sin(th)]; N = nerveBoundaryRadius(q,ray).*ray; % Nerve outline.
sep = r+r'+2*tP+0.05e-3; sep(1:n+1:end) = 0; need = r+tP+0.06e-3; % Extra room (50 and 60 um) so that irregular outlines of the same areas also fit (the shape generator falls back to rounder shapes when needed). % Minimum centre distances and distance to the nerve surface.
for restart = 1:40
    C = (2*rand(rs,n,2)-1).*g.nerveSemiAxes_m*0.7; % Random start.
    for it = 1:2500
        moved = 0; D = C-permute(C,[3 2 1]); % Pairwise differences (n x 2 x n).
        for i = 1:n
            for j = i+1:n
                d = norm(C(i,:)-C(j,:)); if d<sep(i,j), u = (C(i,:)-C(j,:))/max(d,1e-12); s = 0.5*(sep(i,j)-d)+1e-9; C(i,:) = C(i,:)+s*u; C(j,:) = C(j,:)-s*u; moved = max(moved,s); end
            end
            [db,cp,ins] = nearestOnOutline(C(i,:),N); if ~ins, db = -db; end % Signed distance to the nerve surface (positive inside).
            if db<need(i), u = (C(i,:)-cp); u = u/max(norm(u),1e-12); if ~ins, u = -u; end, s = need(i)-db+1e-9; C(i,:) = C(i,:)+s*u; moved = max(moved,s); end
        end
        if moved<1e-9, break; end
    end
    ok = moved<1e-9; for i = 1:n, [db,~,ins] = nearestOnOutline(C(i,:),N); ok = ok && ins && db>=need(i)-1e-9; end % Verify.
    if ok, F = [C r r zeros(n,1)]; return; end
end
error('Could not pack the fascicles inside the nerve.');
end
function [d,cp,ins] = nearestOnOutline(p,P) % Distance, closest point and inside flag of a point with respect to a closed outline.
Q = circshift(P,-1,1); S = Q-P; t = ((p(1)-P(:,1)).*S(:,1)+(p(2)-P(:,2)).*S(:,2))./sum(S.^2,2); t = min(1,max(0,t));
X = P(:,1)+t.*S(:,1); Y = P(:,2)+t.*S(:,2); dist = hypot(p(1)-X,p(2)-Y); [d,i] = min(dist); cp = [X(i) Y(i)]; ins = inpolygon(p(1),p(2),P(:,1),P(:,2));
end
