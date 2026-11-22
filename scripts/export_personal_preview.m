function export_personal_preview() % Make the personal interactive anatomy available before the electrical solve.
[p,a,s]=personalDevice3D(pw6_3d_params(),quicklook3DParams()); % Rebuild the actual planned contacts and patient.
o=struct('params',p,'patient',a,'sites',s,'solvedIds',[],'results',table(),'created',"preview"); % Mark every contact as unsolved.
step=0.1e-3; xs=(-ceil(a.nerveRadius_m/step):ceil(a.nerveRadius_m/step))*step; % Set the transverse display bounds.
model=viewerGeometry(o,xs,xs,(-63:63)*0.2e-3,step,0.2e-3); % Use the same viewer export as the full solve.
model.source="Geometry preview from quicklook3DParams"; % Identify the source separately from the eventual field results.
root=fileparts(fileparts(mfilename('fullpath'))); % Locate this checkout.
folder=fullfile(root,'results','viewer','Cuff60Preview'); % Keep the geometry preview available across future solves.
if ~isfolder(folder), mkdir(folder); end % Prepare the local viewer dataset.
fid=fopen(fullfile(folder,'geometry.js'),'w'); % Open the metadata file.
assert(fid>=0,'Cannot write personal viewer preview.'); % Report output errors.
cleanup=onCleanup(@()fclose(fid)); % Close the dataset after writing.
fprintf(fid,'NV.setModel(%s);\n',jsonencode(model)); % Export geometry without any synthetic field.
end % Return after the interactive preview is ready.
