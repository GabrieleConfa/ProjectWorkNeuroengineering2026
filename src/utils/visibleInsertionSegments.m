function [segments_m,straight_m] = visibleInsertionSegments(strip) % Return the actual continuous single-entry carrier.
segments_m = strip.path_m; % Use the entire collision-checked physical lead.
along = segments_m(:,1:2)*strip.axisXY'; % Measure distance along the nominal insertion direction.
straight_m = [along*strip.axisXY,segments_m(:,3)]; % Keep the same depth for the undeviated comparison.
if isfield(strip,'entryAngle_deg') % Use the actual incoming direction for elastic strips.
    lateral=[-strip.axisXY(2),strip.axisXY(1)]; % Express the local transverse direction in anatomical coordinates.
    straight_m(:,1:2)=straight_m(:,1:2)+(along-along(1))*tand(strip.entryAngle_deg)*lateral; % Anchor the reference at the same entry.
end % Older strips retain their nominal straight reference.
end % Return the two paths in metres.
