function refresh_contact_profiles(matFile) % Render the original potential, axial profile and activation plots with updated personal geometry.
S=load(matFile,'out'); o=S.out; % Read the verified full model.
q=o.params.personalGeometry; f=o.patient.fascicles; % Read the solved patient and fascicle centres.
z=linspace(-q.domain.length_m/2,q.domain.length_m/2,2801)'; % Match the longitudinal sampling of the shared figures.
axesV=zeros(numel(z),height(f),numel(o.solvedIds)); % Allocate axial profiles across all contacts.
T=double(o.mesh.tets); X=o.mesh.nodes_m; % Reuse the same physical mesh.
for k=1:height(f) % Locate each fibre once across all stimulation configurations.
    axesV(:,k,:)=reshape(sampleAxialFields(X,T,o.V_per_A,f.x_m(k),f.y_m(k),z)*1e-6,numel(z),1,[]); % Scale to one microampere.
end % Finish the longitudinal fields.
tr=triangulation(T',X'); geom.footprints=contactFootprints(o.sites); % Recreate the visual footprint metadata from solved pads.
ang=linspace(0,2*pi,300)'; geom.nerveOutline=nerveBoundaryRadius(q,[cos(ang),sin(ang)])'; % Draw the supplied native section.
for sec=unique(o.sites.section)' % Export one sheet per station.
    rows=o.sites.section==sec; cols=find(ismember(o.solvedIds,o.sites.id(rows))); % Locate the six solved contacts.
    target=fullfile(fileparts(matFile),sprintf('Cuff60_3d_section%+d.png',sec)); % Refresh the standard field images.
    plotContactResults(tr,o.V_per_A(:,cols),o.sites(rows,:),o.patient,geom,o.params,o.results(cols,:),target,axesV(:,:,cols)); % Use the updated personal overlays.
    close(gcf); % Release the exported figure before rendering the next station.
end % Finish all original graphical quantities.
end % Return after exporting all profile sheets.
