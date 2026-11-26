function [p, patient, sites] = personalDevice3D(p,q) % Transfer the editable preview geometry to the full conduction model.
%% 1. Use the same provisional anatomy in the mesh and in the contact planner
validateQuicklookAnatomy(q); % Reject overlapping anatomical regions before planning.
f = q.geometry.fascicles_m; % Read the ten editable fascicle outlines in metres.
p.domain.nerveLength_m = q.domain.length_m; % Match the full 60 mm preview.
p.patient.id = "Cuff60"; % Keep personal results separate from the colleague's patient.
R = max(q.geometry.nerveSemiAxes_m); % Read the current outer radius.
F = table(f(:,1),f(:,2),sqrt(f(:,3).*f(:,4)),repmat(q.geometry.tPeri_m,size(f,1),1),'VariableNames',{'x_m','y_m','rEndo_m','tPeri_m'}); % Populate the tissue geometry used by the FEM.
F.a_m=f(:,3); F.b_m=f(:,4); F.angle_deg=f(:,5); % Preserve the supplied oriented ellipses.
nInner = 0; % Small or collinear groups have no strictly interior centres.
if height(F)>2 && rank(f(:,1:2)-mean(f(:,1:2)))==2, nInner=height(F)-numel(unique(convhull(F.x_m,F.y_m))); end % Count centres inside the convex hull for existing plots.
patient = struct('id',p.patient.id,'level',"provisional upper arm",'seed',NaN,'nerveRadius_m',R,'fascicles',F,'generatorInfo',struct('nInner',nInner,'targetInner',NaN,'restarts',0)); % Preserve the existing export format.
sites = loadActiveSites(q.electrode.sitesFile,R); % Load all original measurement rows.
sites.thickness_m(:) = p.electrode.thicknessOverride_m; % Apply the shared finite contact thickness.
%% 2. Transfer the rigid face and the routed internal pad corners
for section = unique(sites.section)' % Build the transverse device at each station.
    cuff = sites.section==section & sites.type=="cuff"; % Identify the two exterior pads.
    rows = find(sites.section==section & sites.type=="intraneural"); % Keep the source ordering of the four internal pads.
    strip = planQuicklookStrip(q,sites(rows,:),sites(cuff,:),"deflect"); % Select a collision-free geometric scenario.
    for k = 1:numel(rows) % Map every active face to its corresponding contact volume.
        face = strip.padFaces_m(:,:,k); % Read the same four vertices shown in the atlas.
        w = face(2,:)-face(1,:); % Compute the in-plane transverse edge.
        w = w/norm(w); % Normalize the pad orientation.
        n = [w(2),-w(1),0]; % Complete the right-handed local coordinate frame.
        i = rows(k); % Find this pad in the FEM table.
        sites.centre_m(i,:) = mean(face,1); % Centre the effective conductive volume on the ideal surface.
        sites.w(i,:) = w; % Store the width direction for mesh and viewer.
        sites.n(i,:) = n; % Store the contact normal for mesh and viewer.
        sites.phi_deg(i) = atan2d(n(2),n(1)); % Rotate the contact box to match the planner.
        sites.r_m(i) = norm(sites.centre_m(i,1:2)); % Update the displaced radial coordinate.
        sites.theta_deg(i) = atan2d(sites.centre_m(i,2),sites.centre_m(i,1)); % Update the displaced polar angle.
    end % Finish the internal contacts at this station.
end % Finish transferring all eight stations.
[~,p.electrode.cuffIndent_m] = cuffContactPlane(q,sites(find(cuff,1),:)); % Set the minimum full-face indentation for the common circular cuff geometry.
for i=find(sites.type=="cuff")' % Calculate each local plane on the supplied oval.
    [~,~,~,~,frame]=cuffContactPlane(q,sites(i,:)); % Preserve the polar landmark and align the face with the local surface tangent.
    sites.centre_m(i,1:2)=frame.centre+sites.thickness_m(i)/2*frame.n; % Place the electrode volume outside its active face.
    sites.n(i,:)=[frame.n 0]; sites.w(i,:)=[frame.t 0]; % Transfer the local orientation to the mesh and viewer.
    sites.phi_deg(i)=atan2d(frame.n(2),frame.n(1)); sites.r_m(i)=norm(sites.centre_m(i,1:2)); % Keep theta_deg as the original landmark angle.
end % Finish positioning external pads.
p.solve3d.sections = []; % Solve each of the 48 active contacts sequentially.
p.solve3d.solver = "pcg"; % Use the memory-efficient iterative solver.
p.personalGeometry = q; % Preserve every provisional setting with the solved mesh.
end % Return the complete patient, geometry settings and contacts.
