function fig = plotCuffApposition(model) % Magnify the local nerve deformation under one finite cuff pad.
%% 1. Read one measured contact and its geometric fit
p = model.parameters; % Read the anatomical and device parameters.
allSites = loadActiveSites(p.electrode.sitesFile,p.geometry.radius_m); % Recover the supplied active-site dimensions.
site = allSites(find(allSites.type=="cuff",1),:); % Select one cuff pad for a readable detail.
[plane_m,depth_m,~,~,frame] = cuffContactPlane(p,site); n=frame.n; t=frame.t; % Fit the rigid pad to the original nerve outline.
angle = linspace(-0.24,0.24,701); % Sample a local angular interval around the pad.
directions = [cosd(site.theta_deg+rad2deg(angle(:))),sind(site.theta_deg+rad2deg(angle(:)))]; % Build transverse rays through the contact area (used by the oblique 3D panel).
%% 1b. Pick a second cuff pad on a different part of the oval, where the flattening is largest
cuffSites = allSites(allSites.type=="cuff",:); % Use only the exterior pads.
depthAll = zeros(size(cuffSites,1),1); % Centre indentation of every cuff pad.
for k=1:size(cuffSites,1), [~,depthAll(k)] = cuffContactPlane(p,cuffSites(k,:)); end % Reuse the shared contact fit.
farFromFirst = abs(mod(cuffSites.theta_deg-site.theta_deg+90,180)-90)>1; % Exclude the first pad and its opposite partner (same anatomical zone).
depthAll(~farFromFirst) = -inf; % Never pick the same zone twice.
[~,iB] = max(depthAll); siteB = cuffSites(iB,:); % Choose the pad with the largest indentation (strongest flattening).
%% 2. Allocate equally scaled views of the nerve and the finite electrode
fig=figure(Name="Cuff–nerve contact",Color="w",Theme="light",Position=[50 55 1220 850]); % Render four panels using millimetres on every spatial axis.
layout=tiledlayout(fig,2,2,TileSpacing="compact",Padding="compact"); % Show anatomical context before the two orthogonal sections.
w=1e3*site.sideXY_m; len=1e3*site.sideZ_m; thick=1e3*p.electrode.thicknessOverride_m; % Use the same contact dimensions as the personal viewer.
%% 3. Show which tissue separates the external pad from the closest fascicle
ax = nexttile(layout); % Place the contact in its anatomical section.
hold(ax,'on'); % Keep the tissue outlines and distance segment on the same axes.
[~,iz] = min(abs(model.nerve.z_m(:,1)-site.z_m)); % Read the actual deformed cross-section at this station.
fill(ax,1e3*model.nerve.x_m(iz,:),1e3*model.nerve.y_m(iz,:),[0.97 0.91 0.77]); % Draw the epineurial region in atlas colours.
theta = linspace(0,2*pi,2049)'; % Sample the outer sheath at sub-degree intervals.
gap = inf; % Retain the minimum distance between this finite pad and any outer perineurium.
for j=1:size(p.geometry.fascicles_m,1) % Compare every fascicle at this station.
    f=p.geometry.fascicles_m(j,:); % Read the centre, semiaxes and rotation.
    rot=[cosd(f(5)) sind(f(5)); -sind(f(5)) cosd(f(5))]; % Rotate the elliptical section into anatomical coordinates.
    outer=([cos(theta) sin(theta)].*(f(3:4)+p.geometry.tPeri_m))*rot+f(1:2); % Sample the outside of the perineurium.
    inner=([cos(theta) sin(theta)].*f(3:4))*rot+f(1:2); % Sample the underlying endoneurial region.
    fill(ax,1e3*outer(:,1),1e3*outer(:,2),[0.65 0.40 0.35],EdgeColor="none"); % Show the outer sheath using the atlas colour.
    fill(ax,1e3*inner(:,1),1e3*inner(:,2),[0.96 0.75 0.68],EdgeColor="none"); % Show the endoneurial interior.
    along=max(-site.sideXY_m/2,min(site.sideXY_m/2,outer*t'-frame.offset)); % Project each sheath sample onto the finite pad.
    padPoints=frame.centre+along*t; % Convert the nearest pad points back into the anatomical frame.
    [distance,idx]=min(vecnorm(padPoints-outer,2,2)); % Measure the shortest sampled gap to this fascicle.
    if distance<gap % Retain the closest pair over all fascicles.
        gap=distance; closest=j; segment=[padPoints(idx,:);outer(idx,:)]; % Store its ID, distance and anatomical endpoints.
    end % Finish the closest-pair update.
end % Finish all outer sheaths.
edge=frame.centre+[-1;1]*site.sideXY_m/2*t; % Locate the two edges of the rigid pad.
plot(ax,1e3*edge(:,1),1e3*edge(:,2),Color=[0.57 0.14 0.83],LineWidth=5); % Identify the external electrode in violet.
plot(ax,1e3*segment(:,1),1e3*segment(:,2),'k--',LineWidth=1.5); % Mark the geometric gap across the epineurium.
[~,~,~,~,frameB]=cuffContactPlane(p,siteB); edgeB=frameB.centre+[-1;1]*siteB.sideXY_m/2*frameB.t; % Locate the second pad shown in the lower-right cut.
plot(ax,1e3*edgeB(:,1),1e3*edgeB(:,2),Color=[0.85 0.40 0.05],LineWidth=5); text(ax,0.80e3*frameB.centre(1),0.80e3*frameB.centre(2),sprintf('AS%d',siteB.id),Color=[0.85 0.40 0.05],HorizontalAlignment="center"); % Mark where the second cut is taken.
text(ax,1e3*p.geometry.fascicles_m(closest,1),1e3*p.geometry.fascicles_m(closest,2),"F"+closest,HorizontalAlignment="center"); % Identify the nearest fascicle.
axis(ax,'equal'); xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]'); % Keep the anatomical proportions in millimetres.
title(ax,sprintf('Pad AS%d to F%d sheath | %.3f mm',site.id,closest,1e3*gap)); % Report the distance in this provisional anatomy.
%% 4. Render the finite pad and the underlying epineurial surface in three dimensions
ax=nexttile(layout); hold(ax,'on'); % Allocate the oblique contact view in pad coordinates.
zLocal=linspace(-0.65,0.65,131)'; % Resolve both axial edges and the return to the unloaded outline.
r=cuffDeformedRadius(p,site,site.z_m+zLocal*1e-3,directions); % Sample the same deformed nerve surface as in the viewer.
across=1e3*(r.*(directions*t')'-frame.offset); height=1e3*(r.*(directions*n')'-plane_m); % Rotate the surface into local transverse and normal coordinates.
surf(ax,across,repmat(zLocal,1,size(r,2)),height,FaceColor=[0.91 0.76 0.53],EdgeColor="none"); % Render the outer epineurium at its actual dimensions.
base=[-w -len;w -len;w len;-w len]/2; % Define the four corners of the active rectangular face.
verts=[base zeros(4,1);base thick*ones(4,1)]; % Extrude the active face through the physical pad thickness.
faces=[1 2 3 4;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8]; % Close the four walls and both faces of the rigid electrode.
patch(ax,Vertices=verts,Faces=faces,FaceColor=[0.57 0.14 0.83],EdgeColor=[0.31 0.08 0.45]); % Display the finite pad above its interface.
axis(ax,'equal'); xlim(ax,[-0.32 0.32]); ylim(ax,[-0.38 0.38]); view(ax,32,21); grid(ax,'on'); % Preserve the physical aspect ratio of the contact.
xlabel(ax,'Across pad (mm)'); ylabel(ax,'Along nerve (mm)'); zlabel(ax,'Height (mm)'); % Use a single unit throughout.
title(ax,sprintf('Rigid pad  |  %.2f × %.2f × %.2f mm',w,len,thick)); % Report width, length and thickness.
%% 5. Compare two transverse cuts at true scale, on two different parts of the oval
for k=1:2 % Same conventions for both pads, so the flattening can be compared directly.
    st = site; if k==2, st = siteB; end % First pad (upper left) and second pad (largest indentation).
    [pl,dep,~,~,fr] = cuffContactPlane(p,st); % Fit the rigid face to the native outline of this pad.
    dirs = [cosd(st.theta_deg+rad2deg(angle(:))),sind(st.theta_deg+rad2deg(angle(:)))]; % Rays through the local contact area.
    nat = nerveBoundaryRadius(p,dirs); dfm = cuffDeformedRadius(p,st,st.z_m,dirs)'; % Native and contacted surface at the pad centre.
    x = 1e3*(dfm.*(dirs*fr.t')-fr.offset); y = 1e3*(dfm.*(dirs*fr.n')-pl); % Contacted surface in the pad frame (mm).
    nx = 1e3*(nat.*(dirs*fr.t')-fr.offset); ny = 1e3*(nat.*(dirs*fr.n')-pl); % Native outline in the same frame (mm).
    ax=nexttile(layout); hold(ax,'on'); % Overlay the native contour, contacted surface and solid pad.
    fill(ax,[x;flipud(x)],[y;-0.1*ones(size(y))],[0.96 0.88 0.72],EdgeColor="none"); % Show tissue below the deformed boundary.
    patch(ax,w*[-1 1 1 -1]/2,[0 0 thick thick],[0.57 0.14 0.83],FaceAlpha=0.55,EdgeColor=[0.42 0.12 0.61]); % Show the pad cross-section with its actual thickness.
    line1=plot(ax,x,y,Color=[0.02 0.44 0.55],LineWidth=2); % Trace the actual pad–nerve interface.
    line2=plot(ax,nx,ny,'k--',LineWidth=0.9); % Show the original nerve outline through the transparent pad.
    axis(ax,'equal'); xlim(ax,[-0.33 0.33]); ylim(ax,[-0.055 thick+0.035]); % Keep both axes at the same magnification.
    xlabel(ax,'Across pad (mm)'); ylabel(ax,'Height (mm)'); % Label the transverse cut.
    title(ax,{sprintf('Transverse cut at AS%d  |  true scale',st.id),sprintf('Surface moves inward by %.2f µm at the pad centre',1e6*dep)}); % Report the imposed indentation.
    text(ax,0,thick/2,sprintf('AS%d pad',st.id),HorizontalAlignment="center"); text(ax,0,-0.035,'Nerve tissue',HorizontalAlignment="center"); % Label pad and tissue.
    legend(ax,[line1 line2],{'Contacted nerve','Native outline'},Location="southoutside",Orientation="horizontal"); % Distinguish the two closely spaced contours.
end % Finish both transverse sections.
title(layout,sprintf('Cuff–nerve contact  |  AS%d (top) and AS%d (largest indentation)',site.id,siteB.id)); % Identify both pads.
end % Return all four views of the same interface.
