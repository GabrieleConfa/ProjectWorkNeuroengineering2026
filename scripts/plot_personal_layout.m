function plot_personal_layout(q) % Regenerate the geometry-dependent views of the shared 3D workflow before meshing.
if nargin<1,q=quicklook3DParams();end % Use the current approved geometry by default.
root=fileparts(fileparts(mfilename('fullpath'))); folder=fullfile(root,'results','3d_preview'); % Save alongside the approved gallery figures.
if ~isfolder(folder),mkdir(folder);end % Prepare a fresh checkout before exporting figures.
[p,patient,sites]=personalDevice3D(pw6_3d_params(),q); % Transfer the exact 48 pad poses into the plotting frame.
sections=unique(sites.section); g=personalStation2D(sections(1),q); % Select the same source station as the original AS2 figures.
fig=figure(Visible="off",Theme="light",Color="white",Position=[60 60 1250 600]); % Prepare the original transverse and axial layout.
t=tiledlayout(fig,1,2,TileSpacing="compact"); ax=nexttile(t); hold(ax,'on'); % Draw the first contacted section at true scale.
for sec=sections' % Draw the actual footprint and carrier for every axial station in the end-on projection.
    station=personalStation2D(sec,q); plot(ax,1e3*station.path(:,1),1e3*station.path(:,2),Color=[0.13 0.38 0.92]); % Show the station-specific insertion.
end % Finish the eight independent sections.
angles=linspace(0,2*pi,241); for k=1:height(patient.fascicles),[xx,yy]=fascicleOutline(patient.fascicles,k,angles);plot(ax,1e3*xx,1e3*yy,'k-');end % Show the anatomical context for all projected sites.
ray=[cos(angles(:)),sin(angles(:))]; rim=ray.*nerveBoundaryRadius(q,ray);plot(ax,1e3*rim(:,1),1e3*rim(:,2),'k-'); % Mark the native envelope of the combined projection.
fp=contactFootprints(sites); % Include every distinct pad footprint across all stations.
for k=1:numel(fp) % Label each contact on the section.
    xy=1e3*fp(k).corners_m([1:4 1],:); fill(ax,xy(:,1),xy(:,2),[0.58 0.22 0.80]); % Draw the finite electrode cross-section.

end % Finish the six contacts.
axis(ax,'equal'); xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); title(ax,'End-on projection · eight stations'); % Label the selected plane.
ax=nexttile(t); hold(ax,'on'); R=max(rim(:,2))*1e3; L=q.domain.length_m*1e3; % Read the radial and longitudinal scales.
for side=[-1 1] % Mark both axial boundaries of the nerve and the cuff exterior.
    plot(ax,[-L L]/2,side*R*[1 1],'k-'); % Show the nerve extent.
    plot(ax,1e3*[-1 1]*p.cuff.length_m/2,side*1e3*(patient.nerveRadius_m+p.cuff.wallThickness_m)*[1 1],Color=[0.2 0.35 0.75]); % Show the extent of the cuff.
end % Finish the reference boundaries.
for kind=["intraneural","cuff"] % Distinguish contact families on the axial map.
    rows=sites.type==kind; h=scatter(ax,1e3*sites.z_m(rows),1e3*sites.centre_m(rows,2),32,'filled',DisplayName=kind); % Plot actual y positions at every station.
    if kind=="intraneural", h.Marker='o'; h1=h; else,h.Marker='s';h2=h;end % Give internal and external contacts distinct symbols.
end % Finish both contact families.
xlabel(ax,'Distance along the nerve, z (mm)'); ylabel(ax,'Height in the cross-section, y (mm)'); title(ax,{'Side view of the 60 mm nerve','Each column is one station; overlapping pads share the same height'}); legend(ax,[h1 h2],Location="southoutside"); % Identify the projection and markers.
title(t,'Current personal device · transverse and axial layout'); exportgraphics(fig,fullfile(folder,'device_layout.png'),Resolution=150);close(fig); % Save the layout without starting a 3D solve.
%% 2. Finite support and contact geometry in the original AS2 region
fig=figure(Visible="off",Theme="light",Color="white",Position=[70 70 1150 610]); ax=axes(fig);hold(ax,'on'); % Inspect the same part of the nerve as the original 3D view.
raw=loadActiveSites(q.electrode.sitesFile,q.geometry.radius_m); shell=buildQuicklookNerve(q,raw(raw.type=="cuff",:)); % Use the contacted surface from the viewer.
sel=shell.z_m(:,1)>=-13e-3 & shell.z_m(:,1)<=-7.5e-3; % Crop the original region of interest to its first two insertions.
surf(ax,1e3*shell.x_m(sel,:),1e3*shell.y_m(sel,:),1e3*shell.z_m(sel,:),FaceColor=[0.86 0.69 0.48],FaceAlpha=0.16,EdgeColor="none"); % Expose the internal device through the tissue.
th=linspace(0,2*pi,121); z=1e3*[min(shell.z_m(sel,1));max(shell.z_m(sel,1))]; % Use the same physical axial crop for fascicles.
for k=1:height(patient.fascicles) % Extrude each fascicle at its current position.
    [x,y]=fascicleOutline(patient.fascicles,k,th); surf(ax,repmat(1e3*x(:)',2,1),repmat(1e3*y(:)',2,1),repmat(z,1,numel(x)),FaceColor=[0.84 0.52 0.43],FaceAlpha=0.12,EdgeColor="none"); % Preserve the rotated shape.
end % Finish the tissue outlines.
for sec=sections' % Select the two stations inside the cropped view.
    rows=sites.section==sec; if sites.z_m(find(rows,1))*1e3<min(z) || sites.z_m(find(rows,1))*1e3>max(z),continue;end % Exclude contacts outside this region.
    strip=planQuicklookStrip(q,raw(rows & raw.type=="intraneural",:),raw(rows & raw.type=="cuff",:),"deflect"); % Use the actual routed lead.
    xyz=1e3*strip.path_m;plot3(ax,xyz(:,1),xyz(:,2),xyz(:,3),Color=[0.13 0.38 0.92],LineWidth=2); % Show its connected carrier.
    for i=find(rows)' % Draw each pad with all three physical dimensions.
        foot=contactFootprints(sites(i,:)); xy=1e3*foot.corners_m; zz=sites.z_m(i)*1e3+[-1 1]*sites.sideZ_m(i)*500; % Extrude the correct in-plane box.
        patch(ax,Vertices=[xy repmat(zz(1),4,1);xy repmat(zz(2),4,1)],Faces=[1 2 3 4;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8],FaceColor=[0.57 0.14 0.83],EdgeColor='none'); % Draw the rigid electrode.
    end % Finish this station.
end % Finish the two leads.
axis(ax,'equal');view(ax,36,23);grid(ax,'on');xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');zlabel(ax,'z (mm)');title(ax,'Device geometry · sections −4 and −3'); % Identify the camera and crop.
exportgraphics(fig,fullfile(folder,'contact_geometry_3d.png'),Resolution=160); % Publish the corrected close-up.
%% 3. Locate the AS2 sampling planes and the axial probes before the volume solve
site=sites(sites.id==2,:); centre=site.centre_m*1e3; er=centre(1:2)/norm(centre(1:2)); % Match the radial longitudinal plane in view_model_3d.
xy=1e3*g.outline; patch(ax,xy(:,1),xy(:,2),centre(3)*ones(size(xy,1),1),[0.10 0.59 0.64],FaceAlpha=0.17,EdgeColor=[0.10 0.45 0.5]); % Show the transverse sampling plane.
r=1e3*nerveBoundaryRadius(q,er); extent=r*[-1 1 1 -1]; % Use the length of the native anatomic diameter in the sampling direction.
patch(ax,extent*er(1),extent*er(2),z([1 1 2 2]),[0.34 0.45 0.78],FaceAlpha=0.15,EdgeColor=[0.28 0.39 0.67]); % Draw the second sampling plane through the contact.
scatter3(ax,centre(1),centre(2),centre(3),38,'r','filled');text(ax,centre(1),centre(2),centre(3),'  AS2'); % Mark the actual source position.
title(ax,'AS2 · transverse and longitudinal sampling planes · geometry only'); % Reserve numerical potentials for the new volume mesh.
exportgraphics(fig,fullfile(folder,'as2_sampling_planes.png'),Resolution=160);close(fig); % Publish the two physical cuts.
fig=figure(Visible="off",Theme="light",Color="white",Position=[60 60 760 660]); ax=axes(fig);hold(ax,'on'); % Draw the probe locations of the 3D field and axial-profile figure.
drawPersonalStation(g); % Match the section that will be sampled by the volume solver.
for k=1:size(g.f,1) % Show one centreline probe per fascicle.
    point=g.f(k,1:2)*1e3; plot(ax,point(1),point(2),'k+');text(ax,point(1),point(2)+0.075,"F"+k,HorizontalAlignment="center"); % Identify the axial sampling line by its transverse position.
end % Finish all probe locations.
fp=contactFootprints(g.sites); % Select the six pads belonging to this sampling plane.
for k=1:numel(fp) % Preserve all six possible source contact locations.
    xy=1e3*fp(k).corners_m;fill(ax,xy(:,1),xy(:,2),[0.57 0.14 0.83]);text(ax,mean(xy(:,1)),mean(xy(:,2))-0.1,"AS"+fp(k).siteIds(1),HorizontalAlignment="center"); % Label the actual pads beside their footprints.
end % Finish the stimulation geometry.
axis(ax,'equal');xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');title(ax,'Station −4 · contact positions and axial sampling lines'); % Display locations in anatomical coordinates.
exportgraphics(fig,fullfile(folder,'station_sampling_lines.png'),Resolution=150);close(fig); % Keep these locations available before the 3D solution.

end
