function vol_m3 = exposedVolume(fx, level_Vpm, current_A)
%EXPOSEDVOLUME Volume of each tissue where |E| exceeds a field level, at a given total current.
%   vol_m3 = exposedVolume(fx, level_Vpm, current_A)
%
%   The FEM is linear, so at current I the field is I times the field per ampere: an element exceeds
%   the level E when its |E| per A is at least E / I. The volume is read from the cumulative curves
%   of fieldExposure (elements sorted by decreasing field).
%
%   Inputs
%     fx         struct from fieldExposure
%     level_Vpm  field level [V/m] (scalar)
%     current_A  total pattern current [A], scalar or vector
%
%   Output
%     vol_m3     numel(current_A) x 3 exposed volume [m^3] per tissue (columns as fx.tissue)

arguments
    fx (1,1) struct
    level_Vpm (1,1) double {mustBePositive}
    current_A (:,1) double {mustBePositive}
end

vol_m3 = zeros(numel(current_A), numel(fx.tissue));
for k = 1:numel(fx.tissue)
    e = fx.E_perA{k};   cv = fx.cumVol_m3{k};
    for j = 1:numel(current_A)
        n = find(e >= level_Vpm / current_A(j), 1, 'last');   % e is sorted in decreasing order
        if ~isempty(n), vol_m3(j, k) = cv(n); end
    end
end
end
