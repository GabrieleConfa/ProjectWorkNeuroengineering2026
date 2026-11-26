function patient = loadPatient(file) % Read a patient from a JSON file. Missing optional data are filled with explicit defaults.
%% Purpose
% The file is the single entry point for patient-specific data. Minimum: id, units "m", the nerve outline, tPeri_m and the
% outline of every fascicle. Optional: function of each fascicle (either "function": "AIN" or a full "composition"), motor
% fraction, motor patches, orientation. Everything optional that is missing is marked as unknown, and the model keeps working.
s = jsondecode(fileread(file)); cat = functionCatalog(); n = numel(s.fascicles);
if ~isfield(s,'source'), s.source = 'measured'; end, if ~isfield(s,'seed'), s.seed = NaN; end, if ~isfield(s,'provenance'), s.provenance = ''; end
if ~isfield(s,'orientation'), s.orientation = struct('x','unknown','y','unknown','status','unknown'); end
if ~isfield(s.nerve,'outline'), th = linspace(0,2*pi,361)'; th(end) = []; a = s.nerve.semiAxes_m(:)'; r = deg2rad(s.nerve.rotation_deg); % Ellipse from semi-axes.
    P = [a(1)*cos(th),a(2)*sin(th)]; s.nerve.outline = P*[cos(r) sin(r);-sin(r) cos(r)]; end
fas = repmat(struct('id',0,'outline',[],'area_m2',0,'roundness',NaN,'composition',[],'motorFraction',NaN,'patches',{{}}),n,1);
for k = 1:n
    if iscell(s.fascicles), f = s.fascicles{k}; else, f = s.fascicles(k); end % jsondecode returns a cell array when the fascicles do not all have the same fields.
    fas(k).id = k; if isfield(f,'id'), fas(k).id = f.id; end
    fas(k).outline = double(f.outline); fas(k).area_m2 = polyarea(fas(k).outline(:,1),fas(k).outline(:,2));
    comp = []; % Functional composition: from "composition", from "function", or unknown.
    if isfield(f,'composition') && isstruct(f.composition) && ~isempty(fieldnames(f.composition)), comp = zeros(1,height(cat)); for i = 1:height(cat), if isfield(f.composition,cat.id(i)), comp(i) = f.composition.(cat.id(i)); end, end
    elseif isfield(f,'function') && ~isempty(f.function), comp = double(cat.id==string(f.function))'; assert(sum(comp)==1,'Fascicle %d: unknown function "%s".',k,string(f.function)); end
    if ~isempty(comp) && sum(comp)>0, comp = comp/sum(comp); end % Normalise (hand-edited files may sum to 0.999 or 1.001).
    fas(k).composition = comp; if isfield(f,'motorFraction'), fas(k).motorFraction = f.motorFraction; end
    if isfield(f,'patches') && ~isempty(f.patches) % jsondecode returns a cell, a 3-D array (equal sizes) or one matrix, depending on the file.
        p = f.patches; if isnumeric(p) && ndims(p)==3, p = arrayfun(@(j)reshape(p(j,:,:),[],2),1:size(p,1),'UniformOutput',false); elseif ~iscell(p), p = {p}; end
        fas(k).patches = cellfun(@double,p,'UniformOutput',false);
    end
end
patient = struct('id',string(s.id),'source',s.source,'seed',s.seed,'units','m','provenance',s.provenance,'orientation',s.orientation, ...
    'nerve',struct('outline',double(s.nerve.outline)),'tPeri_m',s.tPeri_m,'fascicles',fas,'functionScenario','measured','functionSource','measured','functionIds',{cellstr(cat.id)'});
if ~isfield(s.nerve,'semiAxes_m'), patient.nerve.semiAxes_m = [NaN NaN]; else, patient.nerve.semiAxes_m = s.nerve.semiAxes_m(:)'; end
if any(arrayfun(@(f)isempty(f.composition),fas)), patient.functionSource = 'unknown'; end % Some or all functions are not known.
end
