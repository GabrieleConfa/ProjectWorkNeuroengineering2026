function [meshPreview, field] = run_local_mesh_stimulation(p, section) % Mesh and solve one independently testable cuff-contact slab.
arguments % Expose the anatomical and station inputs in one place.
    p (1,1) struct = quicklook3DParams() % Use the same editable inputs as the 3D graphics.
    section (1,1) double = -4 % Select one of the tutor's eight axial sections.
end % Finish input validation.
%% 1. Build the same anatomy and contact used by the preview figures
validateQuicklookAnatomy(p); % Reject invalid oval or fascicle positions before meshing.
sites = loadActiveSites(p.electrode.sitesFile,p.geometry.radius_m); % Read the original 48 active-site rows.
cuff = sites(sites.type=="cuff" & sites.section==section,:); % Select the two finite cuff pads at this station.
assert(height(cuff)==2,'The selected station needs exactly two cuff pads.'); % Keep the local electrical setup explicit.
meshPreview = buildQuicklookMesh(p,cuff); % Mesh the true locally indented contact cross-section.
field = solveQuicklookSlab(p,meshPreview,cuff); % Solve the stationary two-pad volume conductor.
%% 2. Export a mesh view that exposes the finite contact footprints
outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','local_slab'); % Keep generated data outside source files.
if ~exist(outDir,'dir'), mkdir(outDir); end % Create the result directory on its first run.
figMesh = figure(Name="Local contact mesh",Color="w",Position=[80 90 900 650]); % Make the solved tetrahedral slab inspectable.
ax = axes(figMesh); % Select one spatial view.
tr = triangulation(double(field.tets'),field.nodes_m'); % Build the tetrahedral adjacency object.
faces = freeBoundary(tr); % Find the exterior triangular faces of the meshed slab.
patch(ax,Faces=faces,Vertices=1e3*field.nodes_m',FaceColor=[0.83 0.72 0.56],FaceAlpha=0.16, ...
    EdgeColor=[0.38 0.38 0.38],EdgeAlpha=0.18); % Show real tetrahedral exterior edges without hiding the contacts.
hold(ax,'on'); % Overlay both surface electrode node sets.
scatter3(ax,1e3*field.nodes_m(1,field.activeNodes),1e3*field.nodes_m(2,field.activeNodes), ...
    1e3*field.nodes_m(3,field.activeNodes),18,[0.82 0.18 0.12],'filled'); % Mark the source pad's actual mesh nodes.
scatter3(ax,1e3*field.nodes_m(1,field.returnNodes),1e3*field.nodes_m(2,field.returnNodes), ...
    1e3*field.nodes_m(3,field.returnNodes),18,[0.08 0.27 0.58],'filled'); % Mark the grounded return pad nodes.
axis(ax,'equal'); grid(ax,'on'); view(ax,35,28); % Preserve millimetre proportions and allow mesh inspection.
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]'); zlabel(ax,'z [mm]'); % Give the physical mesh axes.
title(ax,sprintf('Section %+d | %d tetrahedra | red: +1 µA, blue: return',section,size(field.tets,2))); % Identify the electrical terminals.
exportgraphics(figMesh,fullfile(outDir,sprintf('mesh_section_%+d.png',section)),Resolution=170); % Save the verified mesh view.
savefig(figMesh,fullfile(outDir,sprintf('mesh_section_%+d.fig',section))); % Keep a rotatable MATLAB mesh figure.
%% 3. Plot the solved potential across the slab centre
extent_m = 1.05*max(p.geometry.nerveSemiAxes_m); % Frame the whole nerve around the local contact.
grid_m = linspace(-extent_m,extent_m,121); % Use a display grid finer than the bulk mesh.
[gx,gy] = meshgrid(grid_m); % Form a transverse section through the pad centre.
query_m = [gx(:),gy(:),repmat(cuff.z_m(1),numel(gx),1)]; % Position all sample points at the solved station.
potential_V = interpolateP1(tr,double(field.V_per_uA_V),query_m); % Interpolate the actual tetrahedral solution.
potential_mV = reshape(1e3*potential_V,size(gx)); % Convert the solved potential to millivolts for display.
figField = figure(Name="Local electrical check",Color="w",Position=[100 100 800 700]); % Create a separate field map.
ax = axes(figField); % Select the cross-section plot.
map = imagesc(ax,1e3*grid_m,1e3*grid_m,potential_mV); % Display the interpolated potential on the sampled plane.
set(map,'AlphaData',double(~isnan(potential_mV))); % Leave points outside the mesh transparent.
set(ax,'Color','w'); % Show outside tissue as white rather than as a false zero-voltage region.
set(ax,'YDir','normal'); hold(ax,'on'); axis(ax,'equal'); % Match the anatomical x-y orientation.
colormap(ax,turbo); colorbar(ax); % Give the voltage scale for the one-microampere solve.
theta = linspace(0,2*pi,241); % Sample all internal fascicle boundaries for visual reference.
for k = 1:size(p.geometry.fascicles_m,1) % Overlay each anatomical perineurium.
    f = p.geometry.fascicles_m(k,:); % Read this fascicle's editable axes and rotation.
    a = f(3)+p.geometry.tPeri_m; b = f(4)+p.geometry.tPeri_m; % Include the outer sheath thickness.
    x = f(1)+a*cos(theta)*cosd(f(5))-b*sin(theta)*sind(f(5)); % Rotate its x boundary.
    y = f(2)+a*cos(theta)*sind(f(5))+b*sin(theta)*cosd(f(5)); % Rotate its y boundary.
    plot(ax,1e3*x,1e3*y,'k-',LineWidth=0.6); % Mark the resistive fascicle sheath.
end % Finish plotting all provisional fascicles.
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]'); % Label the transverse plane.
title(ax,sprintf('Local slab | section %+d | V [mV] for 1 µA | not full-nerve stimulation',section)); % Prevent clinical interpretation.
exportgraphics(figField,fullfile(outDir,sprintf('potential_section_%+d.png',section)),Resolution=170); % Save the solved field map.
savefig(figField,fullfile(outDir,sprintf('potential_section_%+d.fig',section))); % Keep the field map editable in MATLAB.
save(fullfile(outDir,sprintf('station_%+d.mat',section)),'field','p','section','-v7.3'); % Keep the numerical solution for later comparison.
fprintf('Local slab section %+d: %d nodes, %d tetrahedra, return %.6f µA.\n', ...
    section,size(field.nodes_m,2),size(field.tets,2),1e6*field.returnCurrent_A); % Report mesh size and charge balance.
end % Return the verified mesh and local electrical field.
