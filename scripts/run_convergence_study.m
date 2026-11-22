function summary = run_convergence_study(scales, padIds) % Mesh-convergence study of the personal 60 mm model.
%% 0. What this study answers
% The electrical problem: a current I = 1 uA leaves one active pad and flows through anisotropic
% endoneurium, a thin perineurium layer (potential jump), epineurium and muscle, to the grounded
% outer boundary. The potential V solves div(sigma grad V) = 0 (quasi-static volume conductor).
% The finite-element solution depends on the tetrahedra size h. Here the SAME geometry is meshed
% with 3 (or more) element sizes; if the quantities of interest stop changing as h shrinks, the
% mesh is fine enough for them. Every length in the mesh controls is multiplied by one factor.
% Nothing in the anatomy, leads or cuff is changed, and results/main3d is not overwritten.
if nargin<1, scales=[1.96 1.4 1.0 0.714]; end % Size factors relative to the production mesh (ratio 1.4 between levels; 0.714 is finer than production).
if nargin<2, padIds=[1 5 17]; end % AS1 = internal pad, AS5 and AS17 = cuff pads (AS17 has the largest indentation).
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(); % Resolve all sources relative to this file.
outDir=fullfile(root,'results','convergence'); if ~isfolder(outDir),mkdir(outDir);end % Folder for figures and tables.
%% 1. Build the same personal case used by start_nerve
previous=getenv('NERVE_CASE'); restore=onCleanup(@()setenv('NERVE_CASE',previous)); setenv('NERVE_CASE','personal'); %#ok<NASGU> % Select the personal geometry.
[p,patient,sites]=personalDevice3D(pw6_3d_params(),quicklook3DParams()); % Same anatomy, pads and cuff as the production run.
p.solve3d.solver="pcg"; % Memory-efficient iterative solver, as in the production run.
base=p.mesh3d; names=["hContact_m","hEndo_m","hNerve_m","hCoarse_m","hMax_m"]; % Size controls that are scaled together.
%% 2. Mesh and solve at every refinement level
nL=numel(scales); nP=numel(padIds); % Number of meshes and of pads.
R=nan(nL,nP); spread=R; V0=R; Ezmax=R; AF=R; cons=R; N=nan(nL,1); nNodes=N; tMesh=N; % Result arrays.
profile=cell(nL,nP); zAxis=cell(nL,nP); % Axial potential along the probe line.
for k=1:nL % Loop from the coarsest to the finest mesh.
    q=p; for n=names, q.mesh3d.(n)=base.(n)*scales(k); end % Scale all element-size controls by the same factor.
    fprintf('Level %d/%d: size factor %.2f\n',k,nL,scales(k)); % Progress report.
    r=solveContacts(q,patient,sites,padIds); % Mesh, assemble and solve the selected pads.
    N(k)=r.nTets; nNodes(k)=r.nNodes; tMesh(k)=r.time_s; % Mesh statistics.
    R(k,:)=r.R_kOhm; spread(k,:)=r.spread_pct; V0(k,:)=r.Vprobe_V; Ezmax(k,:)=r.EzMax_Vpm; AF(k,:)=r.AF_Vpm2; cons(k,:)=r.cons_pct; % Store the metrics.
    profile(k,:)=r.profile; zAxis(k,:)=r.z; % Store the profiles.
    fprintf('   %d tetrahedra, R = %s kOhm\n',N(k),mat2str(round(r.R_kOhm,3))); % Report access resistances.
end
%% 3. Differences, observed order and Richardson extrapolation
metrics=struct('R_access_kOhm',R,'V_probe_V',V0,'peak_Ez_Vpm',Ezmax,'activating_fn_Vpm2',AF); % Quantities compared.
fn=fieldnames(metrics); rows=table(); r=scales(end-2)/scales(end-1); % Refinement ratio between consecutive levels (order is estimated on the three finest).
for i=1:numel(fn) % Loop over quantities.
    for j=1:nP % Loop over pads.
        f=metrics.(fn{i})(:,j); % Values on the meshes, coarse to fine.
        change=100*abs(diff(f))./abs(f(2:end)); % Relative change between consecutive meshes (%).
        f3=f(end-2:end); % Three finest meshes.
        order=NaN; ext=NaN; gci=NaN; % Defaults for non-monotone sequences.
        if nL>=3 && (f3(1)-f3(2))*(f3(2)-f3(3))>0 % Monotone convergence is required for an order estimate.
            order=log((f3(1)-f3(2))/(f3(2)-f3(3)))/log(r); % Observed order of convergence p.
            ext=f3(3)+(f3(3)-f3(2))/(r^order-1); % Richardson extrapolation to h -> 0.
            gci=100*1.25*abs((f3(3)-f3(2))/f3(3))/(r^order-1); % Grid convergence index (safety factor 1.25).
            if order<0.5 || order>2.5, ext=NaN; gci=NaN; end % Linear tetrahedra cannot converge faster than ~2nd order: outside this range the extrapolation is not trustworthy.
        end
        rows=[rows;table(string(fn{i}),padIds(j),f(end),change(end),order,ext,gci,'VariableNames',{'quantity','pad','finest','changeLastStep_pct','observedOrder','extrapolated','GCI_pct'})]; %#ok<AGROW> % One row per quantity and pad.
    end
end
summary=rows; writetable(rows,fullfile(outDir,'convergence_summary.csv')); % Numbers for the report.
levels=table((1:nL)',scales(:),N,nNodes,tMesh,'VariableNames',{'level','sizeFactor','tetrahedra','nodes','seconds'}); writetable(levels,fullfile(outDir,'convergence_levels.csv')); % Mesh sizes.
% Consistency check against the production solution, which was computed on the same finest mesh.
prod=fullfile(root,'results','main3d','Cuff60_3d.mat'); check=nan(1,nP); % Saved result of start_nerve.
if isfile(prod) && any(abs(scales-1)<1e-12)
    S=load(prod,'out'); [~,ix]=ismember(padIds,S.out.results.id); check=S.out.results.R_kOhm(ix)'; % Saved access resistances.
end
save(fullfile(outDir,'convergence.mat'),'scales','padIds','N','R','spread','V0','Ezmax','AF','cons','profile','zAxis','rows','levels','check'); % Keep everything.
%% 4. Figures
plotConvergence(outDir,scales,padIds,N,metrics,profile,zAxis,rows,check); % Summary and profile figures.
disp(levels); disp(rows); % Print the tables.
end
%% Local functions
function r=solveContacts(p,patient,sites,padIds) % Same mesh, assembly and solver as main_pipeline_3d, for selected pads.
mm=1e-3; tic; % Timer for meshing.
[geom,X,T,elemCell]=meshPersonalDevice(p,patient,sites); % Conforming tetrahedra around every finite pad.
elemTissue=geom.cellTissue(elemCell)'; elemSite=geom.cellSite(elemCell)'; % Tissue and contact label of every element.
checkCuffMeshContact(X,T,elemTissue,elemSite,sites,p.personalGeometry); % Fail if a cuff pad does not share its footprint with the epineurium.
r.time_s=toc; r.nTets=size(T,2); r.nNodes=size(X,2); % Mesh size.
sigT=p.materials.sigmaT_Spm; sigL=p.materials.sigmaL_Spm; Gs=sigT(2)/p.perineurium.realThickness_m; % Conductivities and perineurium sheet conductance.
contactNodes=unique(T(:,elemTissue==5)); ifaceNodes=setdiff(intersect(unique(T(:,elemTissue==1)),unique(T(:,elemTissue==3))),contactNodes); % Fascicle surface nodes.
tl=thinLayerInterface(X,T,elemTissue==1,ifaceNodes(:),Gs); nAll=tl.nNodes; % Duplicate nodes: potential jump across the perineurium.
sIdle=p.electrode.insulatorSigma_Spm; % Idle pads are insulating.
tissueSigma=[sigT(1) sigT(1) sigL(1); 0 0 0; sigT(3) sigT(3) sigL(3); sigT(4) sigT(4) sigL(4); sIdle sIdle sIdle]'; sig=tissueSigma(:,elemTissue); % Anisotropic diagonal conductivities.
[K0,vol]=assembleConductionP1(tl.nodes_m,tl.tets,sig); K0=K0+tl.K; % Stiffness matrix of div(sigma grad V).
Rm=patient.nerveRadius_m; grounded=abs(hypot(X(1,:),X(2,:))-Rm*p.domain.outerRadiusFactor)<1e-10; % Grounded cylindrical boundary.
if p.domain.endBC=="ground", grounded=grounded|abs(abs(X(3,:))-p.domain.nerveLength_m/2)<1e-10; end % Grounded end planes.
gnd=find(grounded); gnd=unique([gnd(:);nonzeros(tl.copyOf(gnd))]); free=true(nAll,1); free(gnd)=false; % Nodes with V = 0.
I0=p.stim.unitCurrent_A; nP=numel(padIds); % Unit current.
[r.R_kOhm,r.spread_pct,r.Vprobe_V,r.EzMax_Vpm,r.AF_Vpm2,r.cons_pct]=deal(zeros(1,nP)); r.profile=cell(1,nP); r.z=cell(1,nP); % Outputs.
for j=1:nP % Solve one pad at a time.
    row=find(sites.id==padIds(j)); es=find(elemSite==padIds(j)); Te=tl.tets(:,es); % Elements of the active pad.
    K=K0+assembleConductionP1(tl.nodes_m,tl.tets(:,es),(p.electrode.contactSigma_Spm-sIdle)*[1;1;1]); % The active pad becomes a good conductor.
    F=accumarray(Te(:),reshape(repmat(I0/sum(vol(es))*vol(es)/4,4,1),[],1),[nAll 1]); % Uniform source: integrates to exactly I0.
    u=zeros(nAll,1); u(free)=solvePCG(K(free,free),F(free),p.solve3d); % Solve K u = F with V = 0 on the outer boundary.
    cn=unique(Te(:)); Vc=mean(u(cn)); % Contact potential.
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
[x,flag,relres]=pcg(A,b,opt.pcgTol,5000,Lc,Lc'); assert(flag==0,'PCG failed (flag %d, residual %.2e).',flag,relres); % Require convergence.
end
function plotConvergence(outDir,scales,padIds,N,metrics,profile,zAxis,rows,check) % Explanatory figures.
nL=numel(scales); nP=numel(padIds); cols=lines(nL); lab=compose('%d tets (h x %.2f)',N,scales(:)); % Colours and legend labels.
kind=@(id)ternary(id==1,'internal pad','cuff pad'); % Pad type for titles.
%% Figure 1: values versus mesh size, with extrapolation
fig=figure(Visible="off",Color="w",Position=[40 40 1250 800],Theme="light"); tl=tiledlayout(fig,2,2,TileSpacing="compact",Padding="compact"); % Four quantities.
labels={'Access resistance (kΩ)','Potential on probe line at pad plane (V per A)','Peak |E_z| on probe line (V/m per A)','Activating function d²V/dz² (V/m² per A)'}; % Axis labels.
fn=fieldnames(metrics); mk='osd^'; % One marker per pad.
for i=1:4 % One panel per quantity.
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); % Axes.
    for j=1:nP % One line per pad, normalized to its own finest value.
        f=metrics.(fn{i})(:,j); plot(ax,N,100*f/f(end),['-' mk(j)],LineWidth=1.6,DisplayName=sprintf('AS%d (%s)',padIds(j),kind(padIds(j)))); % Percentage of the finest-mesh value.
    end
    set(ax,XScale='log'); yline(ax,100,'k:','HandleVisibility','off'); % Reference: finest mesh = 100 %.
    if i==3, ylim(ax,[50 200]); text(ax,0.02,0.95,'values above 200 % clipped (coarse-mesh spike)','Units','normalized','VerticalAlignment','top'); end % Peak field is spiky on coarse P1 meshes.
    xlabel(ax,'Number of tetrahedra'); ylabel(ax,'Value as % of finest mesh'); title(ax,labels{i}); % Labels.
    if i==1, legend(ax,Location='best'); end % One legend.
end
title(tl,'Mesh convergence: each curve should flatten near 100 % as the mesh is refined (production mesh = 3rd point)'); exportgraphics(fig,fullfile(outDir,'convergence_values.png'),Resolution=170); close(fig); % Save.
%% Figure 2: axial profiles on the three meshes
fig=figure(Visible="off",Color="w",Position=[40 40 1250 850],Theme="light"); tl=tiledlayout(fig,nP,2,TileSpacing="compact",Padding="compact"); % Pads x (V, Ez).
for j=1:nP
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); % Potential panel.
    for k=1:nL, plot(ax,zAxis{k,j},profile{k,j},Color=cols(k,:),LineWidth=1.5,DisplayName=lab{k}); end % One curve per mesh.
    xlabel(ax,'z from pad centre (mm)'); ylabel(ax,'V per A (V/A)'); title(ax,sprintf('AS%d (%s): potential on probe line',padIds(j),kind(padIds(j)))); if j==1,legend(ax,Location='best'); end % Labels.
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); % Axial field panel.
    for k=1:nL, plot(ax,zAxis{k,j},-gradient(profile{k,j},zAxis{k,j}*1e-3),Color=cols(k,:),LineWidth=1.5); end % E_z = -dV/dz.
    xlabel(ax,'z from pad centre (mm)'); ylabel(ax,'E_z per A (V/m per A)'); title(ax,sprintf('AS%d: axial electric field',padIds(j))); % Labels.
end
title(tl,'Axial profiles 0.3 mm from the active face, on each mesh'); exportgraphics(fig,fullfile(outDir,'convergence_profiles.png'),Resolution=170); close(fig); % Save.
%% Figure 3: change in the last refinement step
fig=figure(Visible="off",Color="w",Position=[40 40 1000 520],Theme="light"); ax=axes(fig); % Bar chart.
names=["R access","V probe","peak Ez","Act. fn"]; C=zeros(numel(fn),nP); % Matrix quantities x pads.
for i=1:numel(fn), for j=1:nP, C(i,j)=rows.changeLastStep_pct(rows.quantity==string(fn{i}) & rows.pad==padIds(j)); end, end % Last-step changes.
b=bar(ax,C); for j=1:nP, b(j).DisplayName=sprintf('AS%d (%s)',padIds(j),kind(padIds(j))); end % One bar group per quantity.
set(ax,XTickLabel=names); ylabel(ax,'Change between the two finest meshes (%)'); grid(ax,'on'); legend(ax,Location='northwest'); % Labels.
title(ax,'How much each quantity still changes in the last refinement step'); exportgraphics(fig,fullfile(outDir,'convergence_change.png'),Resolution=170); close(fig); % Save.
if all(isfinite(check)), fprintf('Production check (saved results vs this study, scale 1 mesh): max difference %.3g %%\n',max(100*abs(metrics.R_access_kOhm(abs(scales-1)<1e-12,:)-check)./check)); end % Compare with start_nerve output.
end
function v=ternary(c,a,b) % Small inline conditional.
if c, v=a; else, v=b; end
end
