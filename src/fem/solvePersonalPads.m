function r=solvePersonalPads(p,patient,sites,padIds) % Mesh the personal model and solve a chosen set of pads (same method as main_pipeline_3d).
%% Purpose
% Shared by the convergence and edge-refinement studies. It does NOT write to results/main3d.
% Steps: 1) Gmsh mesh, 2) perineurium thin layer, 3) conduction matrix, 4) one PCG solve per pad,
% 5) access resistance and field samples on a probe line 0.3 mm from the active face.
mm=1e-3; tic; % Timer for meshing.
[geom,X,T,elemCell]=meshPersonalDevice(p,patient,sites); % Conforming tetrahedra around every finite pad.
elemTissue=geom.cellTissue(elemCell)'; elemSite=geom.cellSite(elemCell)'; % Tissue and contact label of every element.
checkCuffMeshContact(X,T,elemTissue,elemSite,sites,p.personalGeometry); % Fail if a cuff pad does not share its footprint with the epineurium.
r.time_s=toc; r.nTets=size(T,2); r.nNodes=size(X,2); % Mesh size.
sigT=p.materials.sigmaT_Spm; sigL=p.materials.sigmaL_Spm; Gs=sigT(2)/p.perineurium.realThickness_m; % Conductivities and perineurium sheet conductance.
contactNodes=unique(T(:,elemTissue==5)); ifaceNodes=setdiff(intersect(unique(T(:,elemTissue==1)),unique(T(:,elemTissue==3))),contactNodes); % Fascicle surface nodes.
tl=thinLayerInterface(X,T,elemTissue==1,ifaceNodes(:),Gs); nAll=tl.nNodes; % Duplicate nodes: potential jump across the perineurium.
sIdle=p.electrode.insulatorSigma_Spm; % Idle pads are insulating.
sPDMS=p.electrode.pdmsSigma_Spm; % Conductivity of the PDMS lead bodies and sleeve (insulating by default).
tissueSigma=[sigT(1) sigT(1) sigL(1); 0 0 0; sigT(3) sigT(3) sigL(3); sigT(4) sigT(4) sigL(4); sIdle sIdle sIdle; sPDMS sPDMS sPDMS; sPDMS sPDMS sPDMS]'; sig=tissueSigma(:,elemTissue); % Anisotropic diagonal conductivities.
[K0,vol]=assembleConductionP1(tl.nodes_m,tl.tets,sig); K0=K0+tl.K; % Stiffness matrix of div(sigma grad V).
Rm=patient.nerveRadius_m; grounded=abs(hypot(X(1,:),X(2,:))-Rm*p.domain.outerRadiusFactor)<1e-10; % Grounded cylindrical boundary.
if p.domain.endBC=="ground", grounded=grounded|abs(abs(X(3,:))-p.domain.nerveLength_m/2)<1e-10; end % Grounded end planes.
gnd=find(grounded); gnd=unique([gnd(:);nonzeros(tl.copyOf(gnd))]); free=true(nAll,1); free(gnd)=false; % Nodes with V = 0.
I0=p.stim.unitCurrent_A; nP=numel(padIds); % Unit current.
[r.R_kOhm,r.spread_pct,r.Vprobe_V,r.EzMax_Vpm,r.AF_Vpm2,r.cons_pct,r.endPlanes_pct,r.sideMuscle_pct]=deal(zeros(1,nP)); r.profile=cell(1,nP); r.z=cell(1,nP); % Outputs.
for j=1:nP % Solve one pad at a time.
    row=find(sites.id==padIds(j)); es=find(elemSite==padIds(j)); Te=tl.tets(:,es); % Elements of the active pad.
    K=K0+assembleConductionP1(tl.nodes_m,tl.tets(:,es),(p.electrode.contactSigma_Spm-sIdle)*[1;1;1]); % The active pad becomes a good conductor.
    F=accumarray(Te(:),reshape(repmat(I0/sum(vol(es))*vol(es)/4,4,1),[],1),[nAll 1]); % Uniform source: integrates to exactly I0.
    u=zeros(nAll,1); u(free)=solvePCG(K(free,free),F(free),p.solve3d); % Solve K u = F with V = 0 on the outer boundary.
    cn=unique(Te(:)); Vc=mean(u(cn)); % Contact potential.
    Ig=-(K(gnd,:)*u-F(gnd)); isEnd=abs(abs(tl.nodes_m(3,gnd))-p.domain.nerveLength_m/2)<1e-9; % Current leaving through each grounded node, and which of them lie on the two end planes.
    r.endPlanes_pct(j)=100*sum(Ig(isEnd(:)))/I0; r.sideMuscle_pct(j)=100*sum(Ig(~isEnd(:)))/I0; % Share of the current reaching the grounded end planes and the grounded outer muscle surface.
    r.R_kOhm(j)=Vc/I0/1e3; r.spread_pct(j)=100*(max(u(cn))-min(u(cn)))/Vc; r.cons_pct(j)=100*abs(-sum(K(gnd,:)*u-F(gnd))-I0)/I0; % Access resistance, equipotential check, current balance.
    nrm=sites.n(row,1:2); thick=sites.thickness_m(row); % Pad normal and thickness.
    if sites.type(row)=="cuff", probe=sites.centre_m(row,1:2)-nrm*(thick/2+0.3*mm); else, probe=sites.centre_m(row,1:2)+nrm*0.3*mm; end % Probe 0.3 mm from the active face, inside tissue.
    z=sites.z_m(row)+(-4:0.05:4)'*mm; % Axial line through the pad station.
    v=sampleAxialFields(tl.nodes_m,double(tl.tets),u,probe(1),probe(2),z); assert(all(isfinite(v)),'Probe line leaves the mesh.'); % P1 interpolation along the line.
    Ez=-gradient(v,z); % Axial electric field (V/m).
    h=1*mm; vv=sampleAxialFields(tl.nodes_m,double(tl.tets),u,probe(1),probe(2),sites.z_m(row)+[-h 0 h]'); % Three nodes 1 mm apart, as in the viewer.
    r.Vprobe_V(j)=vv(2); r.EzMax_Vpm(j)=max(abs(Ez)); r.AF_Vpm2(j)=(vv(1)-2*vv(2)+vv(3))/h^2; % Potential at the pad plane, peak |Ez|, activating function d2V/dz2.
    r.profile{j}=v; r.z{j}=(z-sites.z_m(row))/mm; % Profile relative to the pad centre, in mm.
end
end
function x=solvePCG(A,b,opt) % Preconditioned conjugate gradients (same settings as the production solve).
Lc=ichol(A,struct('type','ict','droptol',1e-4,'diagcomp',1e-3)); % Incomplete Cholesky preconditioner.
[x,flag,relres]=pcg(A,b,opt.pcgTol,5000,Lc,Lc'); assert(flag==0 || relres<1e-8,'PCG failed (flag %d, residual %.2e).',flag,relres); % Require convergence (a stagnation just above the tolerance is accepted below 1e-8).
end
