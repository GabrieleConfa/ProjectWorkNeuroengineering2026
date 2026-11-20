function fig = plotCuffApposition(model) % Magnify the local nerve deformation under one finite cuff pad.
%% 1. Read one measured contact and its geometric fit
p = model.parameters; % Read the anatomical and device parameters.
allSites = loadActiveSites(p.electrode.sitesFile,p.geometry.radius_m); % Recover the supplied active-site dimensions.
site = allSites(find(allSites.type=="cuff",1),:); % Select one cuff pad for a readable detail.
n = [cosd(site.theta_deg),sind(site.theta_deg)]; % Read the contact's outward normal.
t = [-n(2),n(1)]; % Read the direction across the contact width.
[plane_m,depth_m] = cuffContactPlane(p,site); % Fit the rigid pad to the original nerve outline.
angle = linspace(-0.24,0.24,701); % Sample a local angular interval around the pad.
directions = [cosd(site.theta_deg+rad2deg(angle(:))),sind(site.theta_deg+rad2deg(angle(:)))]; % Build transverse rays through the contact area.
native_m = nerveBoundaryRadius(p,directions); % Read the original nerve surface along each ray.
deformed_m = cuffDeformedRadius(p,site,site.z_m,directions)'; % Read the contacted surface at the pad centre.
lateral_mm = 1e3*deformed_m.*(directions*t'); % Measure the transverse position in the pad frame.
nativeHeight_um = 1e6*(native_m.*(directions*n')-plane_m); % Express original height above the contact plane.
deformedHeight_um = 1e6*(deformed_m.*(directions*n')-plane_m); % Express the contacted height relative to that plane.
%% 2. Show the transverse contact and the axial transition
fig = figure(Name="Local cuff contact detail",Color="w",Position=[100 100 1050 470]); % Create a dedicated enlargement.
layout = tiledlayout(fig,1,2,TileSpacing="compact",Padding="compact"); % Separate transverse and axial views.
ax = nexttile(layout); % Select the view across the contact width.
hold(ax,'on'); % Overlay original surface, deformed surface and rigid pad.
plot(ax,lateral_mm,nativeHeight_um,Color=[0.50 0.50 0.50],LineStyle="--",LineWidth=1.5); % Show the undeformed nerve.
plot(ax,lateral_mm,deformedHeight_um,Color=[0.13 0.45 0.68],LineWidth=2.5); % Show the nerve touching the pad.
plot(ax,1e3*site.sideXY_m/2*[-1,1],[0,0],Color=[0.04 0.10 0.15],LineWidth=5); % Show the finite 0.20 mm contact.
xlabel(ax,'Across pad [mm]'); ylabel(ax,'Height above pad plane [µm]'); % Label the enlarged dimensions.
title(ax,sprintf('Transverse contact | centre indentation %.2f µm',1e6*depth_m)); % Report the calculated chord depth.
legend(ax,{'Original nerve','Deformed nerve','Rigid pad'},Location="northwest"); % Identify the three physical curves.
z_m = site.z_m+linspace(-0.65e-3,0.65e-3,501)'; % Sample before, across and after the axial pad.
profile_m = cuffDeformedRadius(p,site,z_m,n); % Compute the axial surface at the pad centre.
support_m = nerveBoundaryRadius(p,n); % Read the undeformed radius on that normal.
ax = nexttile(layout); % Select the longitudinal profile view.
plot(ax,1e3*(z_m-site.z_m),1e6*(support_m-profile_m),Color=[0.13 0.45 0.68],LineWidth=2.5); % Show local indentation along z.
hold(ax,'on'); % Add the active axial footprint.
plot(ax,1e3*site.sideZ_m/2*[-1,1],1e6*depth_m*[1,1],Color=[0.04 0.10 0.15],LineWidth=5); % Mark the 0.20 mm pad length.
xlabel(ax,'Distance from pad centre along nerve [mm]'); ylabel(ax,'Surface indentation [µm]'); % Label the axial length and depth.
title(ax,'Axial contact and smooth return'); % Describe the local shoulder.
end % Return the local-contact explanation figure.
