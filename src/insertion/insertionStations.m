function [st, stations] = insertionStations(q, nerve_m)
%INSERTIONSTATIONS One station description per shank of the provided layout, for simulateInsertion.
%   [st, stations] = insertionStations(q, nerve_m)
%
%   Purpose
%     simulateInsertion (Marco's lead-insertion model) works in the plane of one station. This builds
%     the station descriptions exactly as Marco's run_insertion_study does: the nominal insertion axis
%     u points from the centre to the outermost shank pad, the pads keep their distances along that
%     axis, the lead enters at the nerve wall on the -u side and is anchored at the wall on the +u
%     side.
%
%   Inputs
%     q        anatomy in the insertion format (insertionGeometry)
%     nerve_m  (N, 2) nerve outline polygon [m]
%
%   Outputs
%     st        cell (nStations, 1) of structs: u, padS (pad positions along u), tipS (far wall),
%               entryS (entry wall, negative), nerve, section, padHalf (half the pad side)
%     stations  (1, nStations) section numbers (-4 ... +4)

arguments
    q (1,1) struct
    nerve_m (:,2) double
end
raw = loadActiveSites(q.electrode.sitesFile, q.geometry.radius_m);
stations = unique(raw.section)';
st = cell(numel(stations), 1);
for iS = 1:numel(stations)
    sec = stations(iS);
    sites = raw(raw.section == sec & raw.type == "intraneural", :);      % the 4 shank pads
    [~, i] = max(sites.r_m);   u = sites.centre_m(i, 1:2) / sites.r_m(i);  % nominal axis
    if ismember(sec, q.electrode.reverseEntrySections), u = -u; end
    st{iS} = struct('u', u, 'padS', sites.centre_m(:, 1:2) * u', 'tipS', nerveBoundaryRadius(q, u), ...
                    'entryS', -nerveBoundaryRadius(q, -u), 'nerve', nerve_m, 'section', sec, ...
                    'padHalf', sites.sideXY_m(1) / 2);
end
end
