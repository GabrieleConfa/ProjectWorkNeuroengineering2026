function savePatient(patient, file) % Write a patient as a human-readable JSON file (see docs/patient_data.md).
%% Inputs
%   patient  patient struct (loadPatient / virtualPatient / patientFromOutlines format, coordinates in m)
%   file     path of the .json file to write (overwritten); loadPatient(file) reads it back
cat =functionCatalog(); s = patient;
for k = 1:numel(s.fascicles) % The composition is written as {group: fraction} so it can be edited by hand.
    c = s.fascicles(k).composition; o = struct(); if ~isempty(c), for i = 1:height(cat), o.(char(cat.id(i))) = round(c(i),8); end, end
    s.fascicles(k).composition = o;
end
fid = fopen(file,'w'); assert(fid>=0,'Cannot write %s.',file); cleanup = onCleanup(@()fclose(fid)); fwrite(fid,jsonencode(s,PrettyPrint=true));
end
