function patient = buildPatient(p) % The patient used by this run: a file (specific patient) or a virtual patient.
%% Purpose
% Choose the patient in config/quicklook3DParams.m: p.patient.name = "S01" loads data/patients/personal/S01.json (a specific patient
% added by dropping a file); p.patient.file = full path of a JSON file; both empty = reference virtual patient V001
% (p.patient.virtualSeed = 1) with the function hypothesis p.patient.functionScenario (default "none": functions unknown).
%% Inputs / outputs
%   p        parameter struct with p.patient.{name, file, virtualSeed, functionScenario, functionSeed} and the anatomy
%            (p.geometry) used by virtualPatient
%   patient  patient struct (coordinates in m, docs/patient_data.md), already checked by validatePatient (errors if not usable)
pp = p.patient; file = string(pp.file);
if strlength(file)==0 && isfield(pp,'name') && strlength(string(pp.name))>0 % A patient chosen by name: data/patients/personal/<name>.json.
    file = string(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','patients','personal',pp.name+".json"));
    assert(isfile(file),'Patient "%s" not found: put %s in data/patients/personal (see docs/patient_data.md).',pp.name,pp.name+".json");
end
if strlength(file)>0, patient = loadPatient(file);
else, patient = virtualPatient(p,pp.virtualSeed,pp.functionScenario,pp.functionSeed); end
validatePatient(patient,p);
end
