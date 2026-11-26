function validatePatient(patient, q) % Check that a patient (virtual or measured) is usable by the model. Errors explain what to fix.
%% Required content (see docs/patient_data.md)
req = ["id","source","units","nerve","tPeri_m","fascicles"]; for f = req, assert(isfield(patient,f),'Patient is missing the field "%s".',f); end
assert(patient.units=="m" || string(patient.units)=="m",'Coordinates must be in metres (units = "m").');
n = numel(patient.fascicles); assert(n>=1,'The patient has no fascicles.'); cat = functionCatalog();
nerve = polyshape(patient.nerve.outline); assert(area(nerve)>0,'The nerve outline must be a closed polygon.');
inner = polybuffer(nerve,-patient.tPeri_m); polys = cell(n,1);
for k = 1:n
    f = patient.fascicles(k); P = f.outline; assert(size(P,2)==2 && size(P,1)>=8,'Fascicle %d: the outline needs at least 8 points [x y] in metres.',k);
    assert(all(isinterior(inner,P(:,1),P(:,2))),'Fascicle %d lies outside the nerve (or closer than one perineurium thickness to its surface).',k);
    if ~isempty(f.composition), assert(numel(f.composition)==height(cat) && all(f.composition>=0) && abs(sum(f.composition)-1)<1e-4,'Fascicle %d: the composition must have %d non-negative values summing to 1.',k,height(cat)); end
    for j = 1:numel(f.patches), assert(all(isinterior(polyshape(P),f.patches{j}(:,1),f.patches{j}(:,2))),'Fascicle %d: a motor patch lies outside the fascicle.',k); end
    polys{k} = polyshape(P);
end
for i = 1:n-1, for j = i+1:n
    assert(~overlaps(polybuffer(polys{i},patient.tPeri_m),polys{j}),'Fascicles %d and %d overlap (including their perineurium).',i,j);
end, end
end
