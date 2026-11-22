function plotContactResults(tr, V, S, patient, geom, p, res, fileName, axisValues)
%PLOTCONTACTRESULTS For each contact of one station: V on the cross-section through the contact,
% V along z on the axis of every fascicle, and its second difference over one internode (activating-
% function preview), all for 1 uA.
nC = height(S);   R = patient.nerveRadius_m;   f = patient.fascicles;   nF = height(f);
uA = 1e-6;   L = p.domain.nerveLength_m;   dz = p.plot3d.afSpacing_m;
figure(Theme="light",Color="w", Position=[20 20 min(1900, 300*nC + 120) 1000], Name="Contacts");
tl = tiledlayout(3, nC, TileSpacing="loose", Padding="loose");
title(tl, sprintf('Patient %s, station %+d (z = %.1f mm): monopolar, 1 uA', patient.id, S.section(1), ...
      S.z_m(1)*1e3), FontWeight="bold");
n = 181;   [xg, yg] = meshgrid(linspace(-1.15, 1.15, n) * R);
z = linspace(-L/2, L/2, 2801)';
cmap = lines(nF);
X=tr.Points'; T=tr.ConnectivityList'; % Recover the actual solved tetrahedra.
planeAxis=linspace(-1.15,1.15,n)*R; % Use the same transverse samples as the map axes.
mask=false(n,n,2); mask(:,:,1)=true; % Locate only the exact station plane.
[ids,bary]=locateRegularGrid(X,T,{planeAxis,planeAxis,S.z_m(1)+[0 1e-4]},mask); % Reuse one grid location for every active pad.
ok=isfinite(ids); nodes=T(:,ids(ok))'; % Identify samples inside the physical domain.
planeValues=nan(n*n,nC); % Store one cross section per solution.
for k=1:nC % Apply all fields to the common anatomical samples.
    u=double(V(:,k))*uA; map=nan(n); map(ok)=sum(bary(ok,:).*u(nodes),2); % Interpolate the actual potentials.
    planeValues(:,k)=reshape(map',[],1); % Convert ndgrid to the meshgrid convention used below.
end % Finish plane interpolation.
if nargin<9 % Reuse the common axial profiles when rendering from a saved solution.
axisValues=nan(numel(z),nF,nC); % Store longitudinal profiles for the activating function.
for i=1:nF % Sample the true centreline of each fascicle.
    axisValues(:,i,:)=reshape(sampleAxialFields(X,T,V,f.x_m(i),f.y_m(i),z)*uA,numel(z),1,nC); % Preserve all contact columns.
end % Finish computing the longitudinal profiles.
end % Finish the optional sampling step.
for j = 1:nC
    % map at the contact plane
    Vg = planeValues(:,j); % Select the field for this active contact.
    nexttile(tl, j);
    imagesc(1e3*xg(1,:), 1e3*yg(:,1), reshape(1e3*Vg, n, n)); set(gca, YDir="normal"); hold on
    colormap(gca, turbo); cb = colorbar; cb.Label.String = 'V [mV]';
    drawSection(patient, geom, S.id(j), S);
    if isfield(p,'personalGeometry') % Show the physical device on each personal FEM section.
        q = p.personalGeometry; % Recover the input used for this solve.
        raw = loadActiveSites(q.electrode.sitesFile,q.geometry.radius_m); % Read undeviated nominal coordinates.
        same = raw.section==S.section(j); % Select only this axial level.
        cuff = raw(same & raw.type=="cuff",:); % Select its two external electrode faces.
        lead = planQuicklookStrip(q,raw(same & raw.type=="intraneural",:),cuff,"deflect"); % Recreate the physical support.
        plot(1e3*lead.path_m(:,1),1e3*lead.path_m(:,2),'Color',[0.13 0.38 0.92]); % Show the routed lead on the solved potential.
        ang = linspace(0,2*pi,361)'; % Sample the deformed transverse section.
        ray = [cos(ang),sin(ang)]; % Build outward unit directions.
        radii = cuffDeformedRadius(q,cuff,S.z_m(j),ray); % Intersect each ray with the two finite cuff faces.
        plot(1e3*radii(:).*ray(:,1),1e3*radii(:).*ray(:,2),'Color',[0.13 0.38 0.92]); % Mark the local cuff apposition.
    end % Finish the personal station overlay.
    xlim(1.15e3*R*[-1 1]); ylim(1.15e3*R*[-1 1]);
    title(sprintf('AS %d (%s), %.1f kOhm', S.id(j), string(S.type(j)), res.R_kOhm(j)), FontSize=9);
    subtitle(res.location(j), FontSize=7);
    % V(z) along the fascicle axes
    Vz = axisValues(:,:,j); % Select axial profiles for this active contact.
    nexttile(tl, nC + j);
    set(gca, ColorOrder=cmap); hold on
    plot(1e3*z, 1e3*Vz, LineWidth=1);
    xline(1e3*S.z_m(j), 'r:');  xlim(1e3*(S.z_m(j) + [-12 12]*1e-3));
    xlabel('z [mm]'); ylabel('V [mV]'); title('V on the fascicle axes', FontSize=9); grid on
    if j == nC, legend("F" + (1:nF), Location="eastoutside", FontSize=7); end
    % second difference over one internode (activating function preview, V/m^2). The solution is
    % for a positive (anodic) current; for a CATHODIC pulse of 1 uA the sign flips, and a positive
    % value means depolarisation (excitation).
    k = round(dz / (z(2) - z(1)));
    afCath = -(Vz(1:end-2*k, :) - 2*Vz(1+k:end-k, :) + Vz(1+2*k:end, :)) / (k*(z(2)-z(1)))^2;
    nexttile(tl, 2*nC + j);
    set(gca, ColorOrder=cmap); hold on
    plot(1e3*z(1+k:end-k), afCath, LineWidth=1);
    xline(1e3*S.z_m(j), 'r:');  xlim(1e3*(S.z_m(j) + [-12 12]*1e-3));
    xlabel('z [mm]'); ylabel('\Delta^2V/\Deltaz^2 [V/m^2], >0 = depol.');
    title(sprintf('cathodic \\Delta^2V, \\Deltaz = %.1f mm', dz*1e3), FontSize=9); grid on
end
exportgraphics(gcf, fileName);
end
