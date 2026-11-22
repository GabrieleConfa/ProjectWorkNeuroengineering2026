function [geom,X,T,elemCell]=meshPersonalDevice(p,patient,sites) % Generate shared rigid-pad faces before fitting the local nerve surface.
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
z=uniquetol(sort(z),1e-12); radii=cuffDeformedRadius(q,cuff,z,rays); % Compute the contacted outline before volume meshing.
xy=nerveBoundaryRadius(q,rays).*rays; % Retain the undeformed outer scale for mesh sizing.
cfg.outline=xy*1000; cfg.sectionZ=z*1000; cfg.rays=rays; cfg.radii=radii*1000; % Share ordered station geometry with the CAD kernel.
cfg.domain=1000*[p.domain.nerveLength_m,p.cuff.length_m,R*p.domain.outerRadiusFactor,R+p.cuff.wallThickness_m]; % Describe the muscle and insulating cuff.
cfg.fascicles=[1000*[F.x_m F.y_m F.a_m F.b_m] F.angle_deg]; % Preserve each oriented fascicle.
cfg.pads=[sites.id sites.centre_m*1000 1000*[sites.thickness_m sites.sideXY_m sites.sideZ_m] sites.phi_deg]; % Use the real centre and orientation of each box.
cfg.sizes=1000*[m.hContact_m m.hEndo_m m.hNerve_m m.hCoarse_m m.hMax_m m.fineHalfLength_m]; % Refine tissue around the active contacts.
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
