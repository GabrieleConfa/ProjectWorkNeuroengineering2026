function out = simulateInsertion(q, ip, polys, st, rs) % One random insertion of a PDMS lead pushed by forceps through a nerve section.
%% Purpose
% The tip advances along the commanded direction (displacement control: the hand imposes the motion). At each contact
% with a fascicle (outline + perineurium) one of three outcomes is chosen:
%   deflect  : the tip slides along the perineurium (grazing contact, or not enough force to do anything else)
%   displace : the fascicle is pushed aside (needs free room and F_displace <= normal force)
%   pierce   : the tip enters the fascicle (normal force >= puncture force)
% Inputs: q quicklook params, ip insertionParams, polys {n}[Nx2] fascicle outlines (m), st station struct, rs RandStream.
%% 1. Sample the insertion
tPeri = q.geometry.tPeri_m; Rtip = ip.tipDiameter_m/2; Reff = Rtip+tPeri; ds = ip.stepLength_m; % Contact distance = tip radius + perineurium.
n = numel(polys); P = polys; shift = zeros(n,2); % Current outlines and accumulated displacement of every fascicle.
u = st.u(:)'; v = [-u(2) u(1)]; % Nominal axis and its normal.
ang = ip.entryAngleRange_deg*(2*rand(rs)-1); c = cosd(ang); s = sind(ang); d0 = [c*u(1)-s*u(2), c*u(2)+s*u(1)]; % Commanded direction (random angle).
p = st.entryS*u+ip.entryOffsetSD_m*randn(rs)*v; % Entry point: nominal wall position plus a random sideways offset.
nerveIn = polybuffer(polyshape(st.nerve),-Reff); % The tip centre stays inside the nerve.
guard = 0; while ~isinterior(nerveIn,p(1),p(2)) && guard<2000, p = p+ds*d0; guard = guard+1; end % Move to the first point inside the tissue.
Fpush = ip.pushForceMedian_N*exp(ip.pushForceSigma*randn(rs)); % Normal force delivered by the hand in this trial.
Fp = nan(n,1); side = sign(randn(rs,n,1)); side(side==0) = 1; % Puncture force per fascicle (drawn at first contact) and preferred sliding side.
% --- Options added on main (2026-11-26), both off by default so that Marco's force model is unchanged ---
% pierceModel = "probability" (tutor, meeting 2026-11-23; user 2026-11-26): when the tip touches fascicle k it pierces
% it with probability P = P_r(r_k) * P_gamma(gamma) (pierceProbability.m): small fascicles very unlikely, large ones
% likely, and less likely the more the lead hits at a grazing angle; otherwise it bends around it (no displacement).
% Drawn ONCE per fascicle, at its first contact, with the angle of that contact (after a deflection the lead slides
% along the fascicle: drawing again at every step would make piercing almost certain). Law and values: to ask the tutor.
probModel = isfield(ip,'pierceModel') && ip.pierceModel == "probability";
rEq = cellfun(@(x) sqrt(polyarea(x(:,1),x(:,2))/pi), polys(:)); decided = nan(n,1); % Equivalent radius sqrt(A/pi), decisions.
contactAngle = nan(n,1); % Angle (deg) between the heading and the perineurium normal at the first contact with each fascicle.
% returnToLine: in free motion the tip steers back to the straight line of the entry (pure pursuit towards a point
% returnLookahead_m ahead on that line), so that the exit stays on the line of the entry (tutor's rule).
steer = isfield(ip,'returnToLine') && ip.returnToLine; lat0 = dot(p,v); % Lateral position of the entry line.
inner = polybuffer(polyshape(st.nerve),-(tPeri+0.035e-3)); % Region allowed for displaced fascicles.
centre = cell2mat(cellfun(@(x)mean(x,1),polys(:),'UniformOutput',false)); rad = cellfun(@(x)max(vecnorm(x-mean(x,1),2,2)),polys(:)); % Bounding circles for speed.
flagDeflect = false(n,1); flagDisplace = flagDeflect; flagPierce = flagDeflect; % Outcomes touched in this trial.
inside = 0; pierceSeg = zeros(0,3); entryIdx = 0; path = p; reached = false; % State: fascicle being crossed, pierced segments, path, far wall reached.
%% 2. Advance the tip step by step
for it = 1:8000 % Hard cap on the path length (80 mm).
    if ~isinterior(nerveIn,p(1),p(2)) && it>1, reached = dot(p,u)>0; break; end % Wall reached (checked at every step): on the far half the lead is anchored there, on the near half the trial is incomplete.
    hd = d0; % Heading of this step: the commanded direction, or (returnToLine, outside fascicles) towards the entry line.
    if steer && inside==0, tgt = (dot(p,u)+ip.returnLookahead_m)*u + lat0*v; hd = (tgt-p)/norm(tgt-p); end
    pt = p+ds*hd; % Free step along the heading.
    if inside>0 % The tip is inside a fascicle: move freely until it leaves.
        p = pt; path(end+1,:) = p; %#ok<AGROW> % Extend the path.
        if ~inpolygon(p(1),p(2),P{inside}(:,1),P{inside}(:,2)), pierceSeg(end+1,:) = [entryIdx size(path,1) inside]; inside = 0; end %#ok<AGROW> % Record the exit.
        continue
    end
    hit = 0; best = inf; cp = [0 0]; % Find the closest fascicle that the tip would touch.
    for k = 1:n
        if norm(pt-centre(k,:)) > rad(k)+Reff+2*ds, continue; end % Too far to touch.
        [dk,cpk,insideK] = closestOnOutline(pt,P{k}); if insideK, dk = -dk; end % Signed distance to the outline (negative inside).
        if dk < Reff && dk < best, best = dk; hit = k; cp = cpk; end % Closest contact.
    end
    if hit == 0, p = pt; path(end+1,:) = p; continue; end %#ok<AGROW> % Free motion.
    k = hit; nrm = pt-cp; if best<0, nrm = -nrm; end % Outward normal of the perineurium at the contact point.
    nrm = nrm/max(norm(nrm),eps); gamma = acosd(max(-1,min(1,-dot(hd,nrm)))); % Angle between the heading and the normal (0 = head-on).
    if isnan(contactAngle(k)), contactAngle(k) = gamma; end % Record the angle of the first contact.
    Fn = Fpush*cosd(gamma); % Normal force component delivered at the contact.
    if probModel % Pierce with probability P_r(r) * P_gamma(gamma), decided once per fascicle; otherwise bend around it.
        if isnan(decided(k)), decided(k) = rand(rs) < pierceProbability(rEq(k), gamma, ip); end
        if decided(k), action = "pierce"; else, action = "slide"; end
    elseif tand(gamma) > ip.frictionCoeff % Grazing contact: the tip slides around the perineurium.
        action = "slide";
    else % Nearly head-on: compare the cost of displacing and of piercing.
        move = -nrm*(Reff-best); newP = P{k}+move; % Fascicle translated away from the tip by the overlap.
        okRoom = all(inpolygon(newP(:,1),newP(:,2),inner.Vertices(:,1),inner.Vertices(:,2))); % Stays inside the nerve.
        for j = [1:k-1 k+1:n], if okRoom && polyDist(newP,P{j}) < 2*tPeri+ip.fascicleGap_m, okRoom = false; end, end % Keeps the gap to its neighbours.
        Fd = ip.fascicleStiffness_Npm*norm(shift(k,:)+move); % Force to hold the fascicle at its new displacement.
        if okRoom && Fd <= Fn, action = "displace"; P{k} = newP; shift(k,:) = shift(k,:)+move; centre(k,:) = centre(k,:)+move;
        else
            if isnan(Fp(k)), Fp(k) = ip.punctureForceMedian_N*exp(ip.punctureForceSigma*randn(rs)); end % Draw the puncture force once per fascicle.
            if Fn >= Fp(k), action = "pierce"; else, action = "slide"; end % Pierce if the force is enough, otherwise slide.
        end
    end
    switch action
        case "displace", flagDisplace(k) = true; p = pt; path(end+1,:) = p; %#ok<AGROW> % The tip follows the pushed fascicle.
        case "pierce", flagPierce(k) = true; inside = k; entryIdx = size(path,1); p = pt; path(end+1,:) = p; %#ok<AGROW> % The tip enters the fascicle.
        case "slide" % Slide along the tangent of the perineurium with the commanded motion projected on it.
            flagDeflect(k) = true; tg = hd-dot(hd,nrm)*nrm; % Tangential component of the motion.
            if norm(tg) < 0.2, tg = side(k)*[-nrm(2) nrm(1)]; end % Head-on and stuck: the lead buckles sideways to a random side.
            tg = tg/norm(tg); p = p+ds*tg; % Move along the surface.
            [dk,cpk,insideK] = closestOnOutline(p,P{k}); if insideK, dk = -dk; end %#ok<ASGLU> % Re-project onto the contact distance.
            nn = (p-cpk)/max(norm(p-cpk),eps); if insideK, nn = -nn; end
            if dk < Reff, p = p+(Reff-dk)*nn; end % Keep the tip outside the perineurium.
            path(end+1,:) = p; %#ok<AGROW>
    end
    if ~isinterior(nerveIn,p(1),p(2)), reached = dot(p,u)>0; break; end % Wall reached: on the far half the lead is anchored there, on the near half the trial is incomplete.
end
if inside>0, pierceSeg(end+1,:) = [entryIdx size(path,1) inside]; end % Close an open pierced segment.
%% 2b. Finish against the wall: the anchored end of the lead stops with its tip touching the nerve boundary
if reached && inside==0
    while size(path,1)>2 % Retract any overshoot: the tip edge must not cross the wall.
        [dw,~,inN] = closestOnOutline(path(end,:),st.nerve); if ~inN, dw = -dw; end
        if dw >= Rtip-2e-6, break; end
        path(end,:) = [];
    end
    p = path(end,:); pierceSeg(:,2) = min(pierceSeg(:,2),size(path,1)); % Keep the pierced segments valid after the retraction.
    for it = 1:600 % Advance along the commanded direction until the tip edge touches the wall (never past it).
        [dw,~,inN] = closestOnOutline(p,st.nerve); if ~inN, dw = -dw; end
        if dw <= Rtip+2e-6, break; end % Touching (within 2 um): done.
        step = min(ds,dw-Rtip); pt = p+step*d0; blocked = false; % Do not run into a fascicle on the way.
        for k = 1:n
            if norm(pt-centre(k,:)) <= rad(k)+Reff+2*ds && (closestOnOutline(pt,P{k}) < Reff || inpolygon(pt(1),pt(2),P{k}(:,1),P{k}(:,2))), blocked = true; end
        end
        if blocked, break; end
        p = pt; path(end+1,:) = p; %#ok<AGROW>
    end
    [dw,~,inN] = closestOnOutline(p,st.nerve); if ~inN, dw = -dw; end
    reached = dw <= Rtip+2e-6; % Reject trials that do not end against the wall.
else
    reached = false; % The lead must not end inside a fascicle.
end
%% 2c. Entry end: the lead starts against the wall too (tip edge touching the sleeve around the nerve, not beyond it)
if reached
    p1 = path(1,:);
    for it = 1:600
        [dw,~,inN] = closestOnOutline(p1,st.nerve); if ~inN, dw = -dw; end
        if dw <= Rtip+2e-6, break; end % Touching: done.
        step = min(ds,dw-Rtip); pt = p1-step*d0; blocked = false; % Move back along the push direction, never past the touching position.
        for k = 1:n
            if norm(pt-centre(k,:)) <= rad(k)+Reff+2*ds && (closestOnOutline(pt,P{k}) < Reff || inpolygon(pt(1),pt(2),P{k}(:,1),P{k}(:,2))), blocked = true; end
        end
        if blocked, break; end
        p1 = pt; path = [p1;path]; entryIdx = entryIdx+1; pierceSeg(:,1:2) = pierceSeg(:,1:2)+1; %#ok<AGROW> % Prepend the point and keep the pierced indices valid.
    end
    [dw,~,inN] = closestOnOutline(p1,st.nerve); if ~inN, dw = -dw; end
    reached = reached && dw <= Rtip+2e-6; % Reject trials whose entry tip does not touch the wall.
end
%% 2d. Both ends of the lead end exactly ON the wall (the inner face of the wrapped device): the centreline is prolonged to the wall
if reached
    for sideEnd = [1 2] % 1 = entry end, 2 = anchored end.
        if sideEnd==1, p0 = path(1,:); dir = path(1,:)-path(2,:); else, p0 = path(end,:); dir = path(end,:)-path(end-1,:); end
        dir = dir/norm(dir); t = rayWall(p0,dir,st.nerve); % Distance along the heading to the nerve boundary.
        if ~isfinite(t) || t > 0.3e-3, reached = false; break; end % Reject a lead that cannot reach the wall.
        if sideEnd==1, path = [p0+t*dir;path]; entryIdx = entryIdx+1; pierceSeg(:,1:2) = pierceSeg(:,1:2)+1; else, path(end+1,:) = p0+t*dir; end %#ok<AGROW>
    end
end
%% 3. Place the pads along the finished lead (distances measured from the anchored far end)
seg = vecnorm(diff(path),2,2); acum = [0;cumsum(seg)]; Ltot = acum(end); % Arc length along the lead.
a = Ltot-(st.tipS-st.padS(:)); % Each pad keeps its nominal distance from the far end, measured along the lead.
padXY = nan(numel(a),2); padTan = nan(numel(a),2); padTissue = zeros(numel(a),1); % Pad centres, directions and the fascicle they end in.
for k = 1:numel(a)
    if a(k) < 0 || a(k) > Ltot, continue; end % Pad outside the inserted part of the lead.
    padXY(k,:) = interp1(acum,path,a(k)); % Centre of the pad on the path.
    j = find(acum>=a(k),1); j = min(max(j,2),numel(acum)); tg = path(j,:)-path(j-1,:); padTan(k,:) = tg/max(norm(tg),eps); % Local tangent.
    for f = 1:n, if inpolygon(padXY(k,1),padXY(k,2),P{f}(:,1),P{f}(:,2)), padTissue(k) = f; end, end % 0 = epineurium, f = fascicle number.
end
completed = reached; % The lead goes from the entry point to the opposite wall (anchored end).
out = struct('contactAngle_deg',contactAngle,'rEq_m',rEq,'path',path,'polys',{P},'shift',shift,'deflect',flagDeflect,'displace',flagDisplace,'pierce',flagPierce,'pierceSeg',pierceSeg, ...
    'padXY',padXY,'padTan',padTan,'padTissue',padTissue,'Fpush',Fpush,'direction',d0,'completed',completed,'length',Ltot); % Results of this trial.
end
function [d,cp,ins] = closestOnOutline(p,Q) % Distance, closest point and inside flag of a point with respect to a closed outline.
Q2 = circshift(Q,-1,1); S = Q2-Q; L2 = sum(S.^2,2); % Segment vectors and squared lengths.
t = ((p(1)-Q(:,1)).*S(:,1)+(p(2)-Q(:,2)).*S(:,2))./L2; t = min(1,max(0,t)); % Projection parameter on every segment.
X = Q(:,1)+t.*S(:,1); Y = Q(:,2)+t.*S(:,2); dist = hypot(p(1)-X,p(2)-Y); % Closest point of every segment.
[d,i] = min(dist); cp = [X(i) Y(i)]; ins = inpolygon(p(1),p(2),Q(:,1),Q(:,2)); % Smallest distance and inside test.
end
function d = polyDist(A,B) % Minimum distance between two closed outlines (vertex-to-segment, both directions).
d = inf; for i = 1:2
    if i==2, [A,B] = deal(B,A); end
    Q2 = circshift(B,-1,1); S = Q2-B; L2 = sum(S.^2,2)'; t = ((A(:,1)-B(:,1)').*S(:,1)'+(A(:,2)-B(:,2)').*S(:,2)')./L2; t = min(1,max(0,t)); % Projection on every segment.
    d = min(d,min(hypot(A(:,1)-(B(:,1)'+t.*S(:,1)'),A(:,2)-(B(:,2)'+t.*S(:,2)')),[],'all')); % Smallest vertex-to-segment distance.
end
end

function t = rayWall(p,d,P) % Smallest positive distance from p along d to the closed outline P (inf if none).
Q = circshift(P,-1,1); S = Q-P; den = d(1)*S(:,2)-d(2)*S(:,1); w = P-p; % Ray-segment intersection.
tt = (w(:,1).*S(:,2)-w(:,2).*S(:,1))./den; uu = (w(:,1)*d(2)-w(:,2)*d(1))./den; % Parameters along the ray and along each wall segment.
ok = abs(den)>1e-18 & tt>=-1e-12 & uu>=-1e-9 & uu<=1+1e-9; t = min([inf;tt(ok)]);
end
