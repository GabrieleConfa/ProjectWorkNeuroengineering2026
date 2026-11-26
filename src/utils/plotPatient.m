function h = plotPatient(ax, patient, scale, showLabels) % Draw a patient section: fascicles coloured by dominant function, motor patches, labels.
%% Purpose
% One drawing for every patient (virtual or measured). Colour = dominant motor group of the fascicle (see functionCatalog);
% dark patches = motor fibres; fascicles without a known function are grey.
if nargin<3, scale = 1e3; end, if nargin<4, showLabels = true; end
cat = functionCatalog(); hold(ax,'on'); N = patient.nerve.outline; fill(ax,scale*N(:,1),scale*N(:,2),[0.97 0.91 0.77],EdgeColor=[0.35 0.29 0.22],LineWidth=1.3); % Epineurium.
used = false(height(cat),1); h = struct('legend',[],'names',strings(0));
for k = 1:numel(patient.fascicles)
    f = patient.fascicles(k); P = f.outline; col = [0.75 0.75 0.75]; lab = "?"; % Unknown function: grey.
    if ~isempty(f.composition) % Colour = dominant MOTOR group (sensory fibres are 87 percent of the axons and would hide it); label = group and motor fraction.
        mot = f.composition.*(cat.kind=="motor")'; [~,i] = max(mot); col = cat.colour(i,:); used(i) = true; lab = sprintf('%s · motor %.0f%%',cat.id(i),100*f.motorFraction);
    end
    outer = polybuffer(polyshape(P),patient.tPeri_m).Vertices; fill(ax,scale*outer(:,1),scale*outer(:,2),0.6*col,EdgeColor="none"); % Perineurium.
    fill(ax,scale*P(:,1),scale*P(:,2),0.45*col+0.55,EdgeColor=0.6*col,LineWidth=1); % Endoneurium coloured by function.
    for j = 1:numel(f.patches), Q = f.patches{j}; fill(ax,scale*Q(:,1),scale*Q(:,2),[0.55 0.05 0.10],EdgeColor="none",FaceAlpha=0.85); end % Motor patches.
    if showLabels, text(ax,scale*mean(P(:,1)),scale*mean(P(:,2)),sprintf('F%d\n%s',f.id,lab),HorizontalAlignment="center",FontSize=7); end
end
axis(ax,'equal'); hs = gobjects(0); names = strings(0);
for i = find(used)', hs(end+1) = patch(ax,NaN,NaN,0.45*cat.colour(i,:)+0.55,EdgeColor=0.6*cat.colour(i,:)); names(end+1) = cat.id(i)+" · "+cat.name(i); end %#ok<AGROW>
hs(end+1) = patch(ax,NaN,NaN,[0.55 0.05 0.10]); names(end+1) = "motor patches"; h.legend = hs; h.names = names; % Handles for a legend.
end
