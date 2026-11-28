function strip = scenarioStrip(p, sites, cuffSites) %#ok<INUSD> % Lead of the insertion scenario in the format of the original lead planner.
%% Purpose
% Figures, the FEM pad poses and the viewer read lead and pads through the same "strip" structure, whatever planner
% made it. This function returns the scenario lead of the station of the given internal pads.
% Inputs: p (quicklook parameters with p.geometry.scenario from insertionScenario), sites (the four intraneural pads of one
%         station: section, z_m, sideXY_m, sideZ_m [m]), cuffSites (unused, kept for the common call signature).
% Output: strip struct (see planQuicklookStrip): path_m (K x 3 [m]), padFaces_m (4 x 3 x 4 [m]), offset_m, the fascicles of
%         the station after the insertion, the pierced (envelopeFascicles) and displaced fascicle numbers.
sc =p.geometry.scenario; sec = sites.section(1); i = find([sc.stations.section]==sec,1); % Select the station.
assert(~isempty(i),'No scenario for station %+d.',sec); st = sc.stations(i); z = sites.z_m(1); % Station and axial position.
path = [st.path, repmat(z,size(st.path,1),1)]; % Lead centreline with its axial coordinate.
faces = zeros(4,3,height(sites)); % Four vertices of every active face.
for k = 1:height(sites)
    edge = st.padXY(k,:)+[-1;1]*sites.sideXY_m(k)/2*st.padTan(k,:); % Ends of the active face along the lead.
    faces(:,:,k) = [edge([1 2 2 1],:), z+[-1;-1;1;1]*sites.sideZ_m(k)/2]; % Extrude it over the axial pad length.
end
v = [-st.u(2) st.u(1)]; lat = (st.path-st.path(1,:))*v'; [~,peak] = max(abs(lat)); % Largest sideways deviation.
strip = struct('path_m',path,'axisXY',st.u,'entryAngle_deg',st.dirAngle_deg,'section',sec,'offset_m',lat(peak), ...
    'envelopeFascicles',find(st.pierce(:))','padFaces_m',faces,'energy_J',NaN,'EI_Nm2',NaN,'fascicles',{st.polys},'displaced',find(st.displace(:))'); % Same fields as before plus the station fascicles.
end
