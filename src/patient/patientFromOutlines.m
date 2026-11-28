function patient = patientFromOutlines(id, nerveOutline_m, fascicleOutlines_m, tPeri_m, functions) % Build a patient from segmented outlines (histology, MRI, ultrasound).
%% Purpose
% The easiest way to use a real section: export the nerve and the fascicle contours (for example from QuPath or ImageJ) as
% x,y lists in metres, call this function, optionally give one function id per fascicle (see functionCatalog), then savePatient.
% nerveOutline_m: Nx2. fascicleOutlines_m: cell array of Nx2. functions (optional): string array, one group id per fascicle ("" = unknown).
% id: patient name (string). tPeri_m: perineurium thickness [m] (default 14 um, Horn 2023). Output: patient struct in the
% format of loadPatient (source 'measured'); write it with savePatient.
if nargin<4, tPeri_m = 14e-6; end, if nargin<5, functions = strings(numel(fascicleOutlines_m),1); end
cat = functionCatalog(); n = numel(fascicleOutlines_m); fas = repmat(struct('id',0,'outline',[],'area_m2',0,'roundness',NaN,'composition',[],'motorFraction',NaN,'patches',{{}}),n,1);
for k = 1:n
    P = fascicleOutlines_m{k}; fas(k).id = k; fas(k).outline = P; fas(k).area_m2 = polyarea(P(:,1),P(:,2));
    if strlength(functions(k))>0, fas(k).composition = double(cat.id==functions(k))'; assert(sum(fas(k).composition)==1,'Unknown function "%s".',functions(k)); end
end
patient = struct('id',string(id),'source','measured','seed',NaN,'units','m','provenance','Built from segmented outlines.','orientation',struct('x','unknown','y','unknown','status','unknown'), ...
    'nerve',struct('outline',nerveOutline_m,'semiAxes_m',[NaN NaN],'rotation_deg',NaN),'tPeri_m',tPeri_m,'fascicles',fas,'functionScenario','measured','functionSource','measured','functionIds',{cellstr(cat.id)'});
if any(arrayfun(@(f)isempty(f.composition),fas)), patient.functionSource = 'unknown'; end
end
