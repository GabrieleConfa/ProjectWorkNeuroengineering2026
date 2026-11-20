function [segments_m, straight_m] = visibleInsertionSegments(strip) % Show candidate insertions without inventing a through-nerve carrier.
%% 1. Locate the innermost provided active face on each side
path_m = strip.path_m; % Read the exploratory route used by the collision calculation.
axisXY = path_m(end,1:2)-path_m(1,1:2); % Recover the nominal transverse direction.
axisXY = axisXY/norm(axisXY); % Make longitudinal projections independent of nerve size.
centres_m = squeeze(mean(strip.padFaces_m(:,:,1:size(strip.padFaces_m,3)),1)); % Read each finite internal pad centre.
siteAlong_m = centres_m(1,:)*axisXY(1)+centres_m(2,:)*axisXY(2); % Locate the pads along the insertion axis.
halfWidths_m = squeeze(vecnorm(strip.padFaces_m(2,1:2,:)-strip.padFaces_m(1,1:2,:),2,2))/2; % Read half of every active face width.
cutoff_m = min(abs(siteAlong_m)-halfWidths_m(:)'); % Stop each drawn track at the innermost active face.
assert(cutoff_m>0, 'Internal active faces cross the axis; insertion topology needs clarification.'); % Reject an unsupported segment drawing.
%% 2. Leave an explicit gap where the device topology is unknown
along_m = path_m(:,1:2)*axisXY'; % Locate every route sample relative to the nerve axis.
visible = abs(along_m)>=cutoff_m; % Keep only sections supported by the two sets of contact positions.
segments_m = path_m; % Copy the exploratory path before hiding the unsupported middle.
segments_m(~visible,:) = NaN; % Break the plot between the two partial insertion candidates.
straight_m = [along_m*axisXY,repmat(path_m(1,3),size(path_m,1),1)]; % Build the undeviated comparison at the same positions.
straight_m(~visible,:) = NaN; % Never draw the straight reference across the unknown middle.
end % Return two separated schematic segments for plotting only.
