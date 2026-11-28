function [geom,X,T,elemCell]=meshPersonalDevice(p,patient,sites) % Generate shared rigid-pad faces before fitting the local nerve surface.
%% Purpose
% 3D mesher of the Gmsh model (default since D-14). MATLAB prepares the geometry in mm (nerve outline flattened by the cuff
% pads at every transverse profile, irregular fascicles displaced by the leads at each station, PDMS leads and sleeve, the 48
% contact boxes, element sizes), writes it to a JSON file, and runs src/mesh/mesh_personal.py (Python + gmsh/OpenCASCADE,
% interpreter from NERVE_GMSH_PYTHON or "python") through system(). The Python script builds the solids, fragments them
% into conforming regions, meshes with linear tetrahedra and writes nodes, elements and the region of each element back.
%% Inputs / outputs
%   p         project parameters with p.mesh3d (element sizes hContact_m ... hMax_m, fineHalfLength_m [m]), p.domain, p.cuff
%             and p.personalGeometry (the Gmsh-model parameters q with the insertion scenario)
%   patient   patient struct (fascicle table with x_m, y_m, a_m, b_m, angle_deg; nerveRadius_m)
%   sites     contact table (centre_m, thickness_m, sideXY_m, sideZ_m [m], phi_deg)
%   geom      struct describing the regions (cellTissue, cellSite, cellFascicle: tissue code / contact id / fascicle of
%             each region), footprints and planes_m for the plots
%   X         3 x nNodes node coordinates [m];  T  4 x nTets node indices;  elemCell  region of each element
%% 1. Export the device in millimetres for OpenCASCADE
R=patient.nerveRadius_m; F=patient.fascicles; m=p.mesh3d; % Read the user-configurable anatomy and element sizes.
cuff=sites(sites.type=="cuff",:); q=p.personalGeometry; % Read the external contacts at their real axial coordinates.
angles=(0:63)'*2*pi/64; corners=[]; z=[-1;1]*p.domain.nerveLength_m/2; % Include both ends of the 60 mm domain.
for j=1:height(cuff) % Add measured stations, pad edges and smooth recovery shoulders.
    pad=cuff(j,:); [~,~,~,~,frame]=cuffContactPlane(q,pad); % Use the same flat face as the atlas.
    corners=[corners;mod(frame.cornerAngles,2*pi)]; %#ok<AGROW> % Include the two side edges exactly.
    h=pad.sideZ_m/2; blend=q.electrode.cuffTransition_m; % Read the axial pad and shoulder lengths.
    z=[z;pad.z_m+[-h-blend;-h-blend/2;-h;0;h;h+blend/2;h+blend]]; %#ok<AGROW> % Build a transverse profile at each physical transition.
end % Finish collecting matching transverse profiles.
near=any(abs(atan2(sin(angles-corners'),cos(angles-corners')))<0.012,2); % Avoid sliver edges next to an exact pad corner.
angles=uniquetol(sort([angles(~near);corners]),1e-10); rays=[cos(angles) sin(angles)]; % Reuse the same ordered vertex directions on all sections.
sc=[]; if isfield(q.geometry,'scenario') && ~isempty(q.geometry.scenario), sc=q.geometry.scenario; end % Random fascicles, PDMS leads and sleeve (insertion mode).
z=uniquetol(sort(z),1e-12); radii=cuffDeformedRadius(q,cuff,z,rays); % Compute the contacted outline before volume meshing.
xy=nerveBoundaryRadius(q,rays).*rays; % Retain the undeformed outer scale for mesh sizing.
cfg.outline=xy*1000; cfg.sectionZ=z*1000; cfg.rays=rays; cfg.radii=radii*1000; % Share ordered station geometry with the CAD kernel.
cfg.domain=1000*[p.domain.nerveLength_m,p.cuff.length_m,R*p.domain.outerRadiusFactor,R+p.cuff.wallThickness_m]; % Describe the muscle and insulating cuff.
cfg.fascicles=[1000*[F.x_m F.y_m F.a_m F.b_m] F.angle_deg]; % Preserve each oriented fascicle.
cfg.pads=[sites.id sites.centre_m*1000 1000*[sites.thickness_m sites.sideXY_m sites.sideZ_m] sites.phi_deg]; % Use the real centre and orientation of each box.
cfg.sizes=1000*[m.hContact_m m.hEndo_m m.hNerve_m m.hCoarse_m m.hMax_m m.fineHalfLength_m]; % Refine tissue around the active contacts.
if isfield(m,'hEdge_m') && m.hEdge_m>0, cfg.edgeSize=1000*m.hEdge_m; end
if isfield(m,'hEdgePads'), cfg.edgePads=m.hEdgePads(:)'; end % Optional list of pad IDs to refine (default: all pads).
if ~isempty(sc) % Export irregular fascicles, PDMS lead bodies and the oval sleeve for Gmsh (millimetres).
    L=p.domain.nerveLength_m; nF=numel(sc.polys0); cfg.fascZ=cell(1,nF); cfg.fascXY=cell(1,nF); % Per fascicle: axial levels and one closed polygon at each level.
    for k=1:nF % An undisturbed fascicle is a straight prism (two end levels); a displaced one needs its fade-out levels at each station.
        zl=[-L/2;L/2]; for s=sc.stations, if any(s.shift(k,:)~=0), zl=[zl;s.z_m+sc.bump_m*[-1;0;1]]; end, end % Levels where this fascicle moves.
        zl=uniquetol(sort(zl),1e-12); cfg.fascZ{k}=1000*zl(:)'; cfg.fascXY{k}=cell(1,numel(zl)); % Axial levels of fascicle k.
        for i=1:numel(zl) % Outline at every level.
            P=scenarioFascicles(sc,zl(i)); idx=round(linspace(1,size(P{k},1)+1,49)); idx(end)=[]; cfg.fascXY{k}{i}=1000*P{k}(idx,:); % 48 vertices per fascicle in the mesh.
        end
    end
    cfg.leads=cell(1,numel(sc.stations)); % PDMS body of each lead: the buffered centreline, extruded over its axial thickness.
    for i=1:numel(sc.stations)
        pb=polybuffer(sc.stations(i).path,'lines',sc.leadWidth_m/2,'JointType','miter'); pb=intersect(pb,polyshape(sc.nerve)); rg=regions(pb); [~,big]=max(arrayfun(@area,rg)); xy=rg(big).Vertices; xy=xy(all(isfinite(xy),2),:); % Band of the lead clipped flush with the nerve wall: it ends against the wrapped device, not beyond it.
        cfg.leads{i}=struct('z',1000*sc.stations(i).z_m,'thickness',1000*sc.leadThickness_m,'xy',1000*xy); % Band outline and axial extent.
    end
    th=linspace(0,2*pi,97)'; th(end)=[]; ray=[cos(th),sin(th)]; r0=nerveBoundaryRadius(q,ray); pts=r0(:).*ray; % Native oval outline (96 vertices).
    t=circshift(pts,-1,1)-circshift(pts,1,1); nrm=[t(:,2),-t(:,1)]./vecnorm(t,2,2); % Outward normal at every vertex.
    cfg.sleeveXY={1000*(pts+sc.sleeveThickness_m*nrm),1000*(pts+sc.sleeveThickness_m*nrm)}; cfg.sleeveZ=1000*p.cuff.length_m/2*[-1 1]; % Oval PDMS sleeve: a straight prism between the cuff ends (the 4 to 9 um pad indentations are negligible against its thickness).
end
%% 2. Generate the conforming tetrahedra and label every cell
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); % Isolate concurrent runs and remove intermediates.
input=fullfile(folder,'device.json'); output=fullfile(folder,'mesh.mat'); % Keep both files local to this run.
writelines(jsonencode(cfg),input); % Serialize dimensions without losing the pose of the pads.
python=getenv('NERVE_GMSH_PYTHON'); % Allow another installation to provide gmsh, scipy and numpy.
if isempty(python) % Discover a repository-local virtual environment first.
    root=fileparts(fileparts(fileparts(mfilename('fullpath')))); % Locate the project on any operating system.
    if ispc, executable=fullfile('Scripts','python.exe'); fallback='python';else, executable=fullfile('bin','python');fallback='python3';end % Match virtual-environment layouts.
    choices={fullfile(root,'.venv',executable),fullfile(getenv('HOME'),'.cache','nerve-gmsh',executable)}; % Also support the existing Mac installation.
    python=fallback; for k=1:numel(choices),if isfile(choices{k}),python=choices{k};break;end,end % Fall back to the Python on PATH.
end % An explicit NERVE_GMSH_PYTHON always takes precedence.
script=fullfile(fileparts(mfilename('fullpath')),'mesh_personal.py'); % Locate the matching meshing script.
status=system(sprintf('"%s" "%s" "%s" "%s"',python,script,input,output)); % Generate and classify the tetrahedra.
assert(status==0,'Gmsh could not generate the personal device mesh.'); % Reject a partial export.
data=load(output); X=data.X; T=double(data.T); elemCell=double(data.elemCell); % Adapt the connectivity to the shared MATLAB solver.
geom=struct('cellTissue',double(data.roles(:,1)),'cellFascicle',double(data.roles(:,2)),'cellSite',double(data.roles(:,3))); % Keep the solver's existing region identifiers.
geom.cellZone=zeros(size(geom.cellSite)); geom.sectionZ_m=z; % Store the reference edges for local surface recovery.
geom.footprints=contactFootprints(sites); geom.planes_m=unique(sites.z_m)'; % Supply the physical pad positions to the section plots.
th=linspace(0,2*pi,300)'; geom.nerveOutline=nerveBoundaryRadius(p.personalGeometry,[cos(th) sin(th)])'; % Preserve the supplied native outline.
end
