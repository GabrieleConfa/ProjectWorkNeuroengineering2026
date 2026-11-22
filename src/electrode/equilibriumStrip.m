function strip=equilibriumStrip(p,sites,cuff) % Equilibrium of a planar strip in a fixed anatomical section.
%% 1. Geometry, bending stiffness and unloaded insertion imperfection
e=p.electrode; f=p.geometry.fascicles_m; % Read the current strip and fascicle dimensions.
assert(all([e.stripWidth_m,e.stripThickness_m,e.stripYoungModulus_Pa,e.unsupportedLength_m,e.foundation_Npm2]>0),'Positive strip dimensions and stiffness required.'); % Reject invalid mechanical inputs.
[~,i]=max(sites.r_m); u=sites.centre_m(i,1:2)/sites.r_m(i); % Recover the nominal pad axis from the device coordinates.
if isfield(e,'reverseEntrySections') && ismember(sites.section(1),e.reverseEntrySections), u=-u; end % Select the entry end without moving the nominal pad coordinates.
v=[-u(2),u(1)]; % Orient the lateral deflection axis relative to insertion.
padS=sites.centre_m(:,1:2)*u'; tip=max(padS+sites.sideXY_m/2)+e.tipOverhang_m; % Keep the distal tip beyond all four pads.
entry=-cuffDeformedRadius(p,cuff,sites.z_m(1),-u); % Start exactly at the local nerve surface.
assert(e.tipOverhang_m>0 && tip<nerveBoundaryRadius(p,u),'Invalid lead length.'); % Reject a through-nerve device.
N=101; s=linspace(entry,tip,N)'; dx=s(2)-s(1); z=sites.z_m(1); % Discretize the insertion axis in metres.
EI=e.stripYoungModulus_Pa*e.stripThickness_m*e.stripWidth_m^3/12; % In-plane bending about z uses the in-plane width cubed.
D=diff(eye(N),2)/dx^2; K=EI*dx*(D'*D); % Discretize the small-slope Euler–Bernoulli bending energy.
H=K+e.foundation_Npm2*dx*eye(N); % Add the assumed linear restoring force from surrounding tissue.
stream=RandStream('mt19937ar','Seed',e.routeSeed+100+sites.section(1)); % Give each insertion a reproducible unloaded imperfection.
t=(s-entry)/(tip-entry); y0=e.initialBow_m*t.^2.*(2*rand(stream)-1+(2*rand(stream)-1)*sin(pi*t)); % Perturb the unloaded strip only.
%% 2. Apply the station-specific entry angle about the fixed insertion point
angleRange=0; % Reproduce straight entry for saved parameter sets predating angular variation.
if isfield(e,'entryAngleRange_deg'), angleRange=e.entryAngleRange_deg; end % Use the configured half-width for new runs.
validateattributes(angleRange,{'numeric'},{'scalar','real','finite','nonnegative','<',90}); % Exclude backward entry directions.
angleDeg=angleRange*(2*rand(stream)-1); % Sample the insertion angle symmetrically before evaluating any obstacle.
axisOffset=(s-entry)*tand(angleDeg); % Define the unloaded straight insertion line through the original entry.
load=K*y0+e.foundation_Npm2*dx*axisOffset; % Restore the lead toward its incoming line in the absence of contact.
cx=f(:,1:2)*u'; cy=f(:,1:2)*v'; % Project fascicle centres into the local insertion frame.
r=max(f(:,3:4),[],2)+p.geometry.tPeri_m; % Enclose each rotated perineurium in a conservative disk.
radius=r+e.stripWidth_m/2; % Include the finite width of the carrier.
%% 3. Enumerate feasible passages and solve their constrained equilibria
lo=-e.maxDeflection_m*ones(N,1); hi=-lo; % Set the available lateral travel.
ang=p.geometry.nerveRotation_deg; rot=[cosd(ang) -sind(ang);sind(ang) cosd(ang)]; % Express the insertion frame along anatomical axes.
axes=p.geometry.nerveSemiAxes_m; U=(u*rot)./axes; V=(v*rot)./axes; % Write the elliptical wall as a quadratic in lateral displacement.
mid=-s*dot(U,V)/dot(V,V); half=sqrt(max(0,mid.^2+(1-s.^2*dot(U,U))/dot(V,V))); % Find both intersections at every depth.
lo=max(lo,mid-half); hi=min(hi,mid+half); % Enforce tissue containment during the energy minimization.
lo(1:2)=axisOffset(1:2); hi(1:2)=lo(1:2); % Clamp entry position and initial insertion direction.
best=Inf; strip=[]; % Retain the lowest energy across all admissible passages.
visit(1,lo,hi,false(size(r))); % Compare both sides of every reachable fascicle.
assert(~isempty(strip),'No feasible single-entry lead at station %+d for these inputs.',sites.section(1)); % Report infeasible anatomy explicitly.
    function visit(j,L,U,puncture) % Bound each candidate passage before solving the beam.
        if any(L>U), return; end % This corridor is geometrically closed.
        if j<=numel(r) % Add a perineurial contact constraint.
            ds=max(abs(s-cx(j))-dx,0); near=ds<radius(j); % Bound the full interval between consecutive samples.
            h=sqrt(max(0,radius(j)^2-ds.^2)); % Compute conservative upper and lower exclusion envelopes.
            if ~any(near) || all(L(near)>=cy(j)+h(near)) || all(U(near)<=cy(j)-h(near)) % This obstacle lies outside the available corridor.
                visit(j+1,L,U,puncture); return % Move on without duplicating the same shape.
            end % Finish the reachability check.
            upper=L; upper(near)=max(L(near),cy(j)+h(near)); % Constrain the beam above this fascicle.
            lower=U; lower(near)=min(U(near),cy(j)-h(near)); % Constrain the beam below it.
            visit(j+1,upper,U,puncture); visit(j+1,L,lower,puncture); % Compare both elastic equilibria.
            if isfinite(e.appliedInsertionForce_N) && isfinite(e.calibratedPunctureForce_N) && e.calibratedPunctureForce_N>0 % Require device-specific evidence for rupture.
                puncture(j)=true; visit(j+1,L,U,puncture); % Test a crossing using the actual incidence and available load below.
            end % Finish conditional rupture.
            return % All descendant passages have now been evaluated.
        end % Every fascicle has an assigned admissible passage.
        y=boundedBeam(H,load,L,U); % Solve stationary elastic equilibrium in this corridor.
        energy=0.5*EI*dx*sum((D*(y-y0)).^2)+0.5*e.foundation_Npm2*dx*sum((y-axisOffset).^2); % Compare total elastic energies in joules.
        if energy>=best, return; end % A positive puncture work can only increase this cost.
        xy=s*u+y*v; tang=gradient(y,dx); % Recover the actual centreline and slopes.
        rr=vecnorm(xy,2,2); boundary=cuffDeformedRadius(p,cuff,z,xy./max(rr,eps)); % Check the deformed nerve at every point.
        if any(rr>boundary(:)+1e-10), return; end % Reject paths leaving the nerve.
        faces=zeros(4,3,height(sites)); % Allocate all finite active pad faces.
        for k=1:height(sites) % Place each pad along the deformed support.
            centre=padS(k)*u+interp1(s,y,padS(k))*v; % Preserve the nominal insertion coordinate.
            tangent=u+interp1(s,tang,padS(k))*v; tangent=tangent/norm(tangent); % Align the face with the local strip.
            edge=centre+[-1;1]*sites.sideXY_m(k)/2*tangent; % Preserve the true transverse active length.
            direction=edge(2,:)-edge(1,:); off=f(:,1:2)-edge(1,:); % Measure distance to the whole pad.
            proj=max(0,min(1,off*direction'/sum(direction.^2))); % Include both finite endpoints.
            if any(vecnorm(off-proj.*direction,2,2)<=r+e.stripThickness_m/2 & ~puncture), return; end % Reject contacts intersecting protected sheaths.
            edgeR=vecnorm(edge,2,2); limit=cuffDeformedRadius(p,cuff,z,edge./max(edgeR,eps)); % Check pad containment.
            if any(edgeR>limit(:)+1e-10), return; end % Reject active faces outside the nerve.
            faces(:,:,k)=[edge([1 2 2 1],:),z+[-1;-1;1;1]*sites.sideZ_m(k)/2]; % Extrude the active face over its axial length.
        end % Finish all pad checks.
        crossed=[]; % Keep actual crossed sheaths for this candidate.
        for k=find(puncture)' % Check actual entry into each potentially penetrable sheath.
            dsxy=diff(xy); rel=xy(1:end-1,:)-f(k,1:2); % Define segment–sheath intersection quadratics.
            aa=sum(dsxy.^2,2); bb=sum(rel.*dsxy,2); cc=sum(rel.^2,2)-r(k)^2; % Use the physical perineurial envelope.
            root=(-bb-sqrt(max(0,bb.^2-aa.*cc)))./aa; % Find the first contact along each segment.
            hit=find(bb.^2>=aa.*cc & root>=0 & root<=1,1); % Locate the first actual sheath crossing.
            if isempty(hit), return; end % Require a real puncture before waiving the exclusion envelope.
            normal=(rel(hit,:)+root(hit)*dsxy(hit,:))/r(k); % Evaluate the local surface normal.
            cosine=abs(dsxy(hit,:)*normal'/sqrt(aa(hit))); % Project the load onto that normal.
            buckling=pi^2*e.stripYoungModulus_Pa*min(e.stripWidth_m*e.stripThickness_m^3,e.stripThickness_m*e.stripWidth_m^3)/(48*e.unsupportedLength_m^2); % Ideal free-column capacity limited by its weaker axis.
            if min(e.appliedInsertionForce_N,buckling)*cosine<e.calibratedPunctureForce_N, return; end % Reject insufficient normal force.
            energy=energy+e.calibratedPunctureForce_N*p.geometry.tPeri_m/max(cosine,eps); % Approximate perforation work across the sheath.
            crossed(end+1)=k; % Record the conditionally penetrated fascicle.
        end % Finish the calibrated-force screen.
        if energy>=best, return; end % Compare total elastic and puncture work.
        best=energy; [~,peak]=max(abs(y)); % Record the chosen equilibrium and signed maximum deviation.
        strip=struct('path_m',[xy,repmat(z,N,1)],'axisXY',u,'entryAngle_deg',angleDeg,'section',sites.section(1),'offset_m',y(peak),'envelopeFascicles',crossed,'padFaces_m',faces,'energy_J',energy,'EI_Nm2',EI); % Return one consistent device for figures and FEM.
    end % Finish passage evaluation.
end % Return the minimum-energy feasible strip.
function y=boundedBeam(H,f,L,U) % Active-set solution of a positive-definite box-constrained quadratic.
y=max(L,min(U,zeros(size(L)))); locked=L==U; % Start from a feasible point and fix the clamped entry nodes.
lower=y==L; upper=y==U; % Track currently constrained lateral positions.
for iter=1:2000 % Release or add constraints until the force residual meets the KKT conditions.
    fixed=lower|upper; free=~fixed; trial=y; % Hold the current contact nodes fixed.
    trial(free)=H(free,free)\(f(free)-H(free,fixed)*y(fixed)); % Solve equilibrium on the unconstrained nodes.
    d=trial-y; ratio=inf(size(y)); % Find the first boundary intercepted on the Newton step.
    down=free & trial<L; up=free & trial>U; % Identify constraints reached along this step.
    ratio(down)=(L(down)-y(down))./d(down); ratio(up)=(U(up)-y(up))./d(up); % Compute feasible step fractions.
    [step,k]=min(ratio); % Select the first blocking node.
    if step<1 % Add its contact reaction before continuing.
        y=y+max(0,step)*d; lower(k)=down(k); upper(k)=up(k); y(k)=max(L(k),min(U(k),trial(k))); continue % Keep the iterate feasible.
    end % The full free step fits the corridor.
    y=trial; grad=H*y-f; wrong=zeros(size(y)); % Calculate residual nodal force.
    wrong(lower & ~locked)=-grad(lower & ~locked); wrong(upper & ~locked)=grad(upper & ~locked); % Test reaction force signs.
    [err,k]=max(wrong); if err<1e-11, return; end % Stop when no contact reaction points the wrong way.
    lower(k)=false; upper(k)=false; % Release the most violated constraint.
end % Failure to reach stationary equilibrium is a solver error.
error('Elastic contact equilibrium failed to converge.'); % Do not export an unconverged trajectory.
end % Finish the constrained beam solver.
