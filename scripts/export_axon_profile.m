function export_axon_profile(matFile,siteId,fascicleId,outputFile) % Export the extracellular drive of one axial fibre for NEURON.
%% 1. Read the solved contact and its anatomical position
S = load(matFile,'out'); % Read the actual 3D FEM solution.
o = S.out; % Access the mesh, geometry and potentials together.
j = find(o.solvedIds==siteId,1); % Locate the requested solution column.
assert(~isempty(j),'Contact was not solved.'); % Reject unavailable fields.
site = o.sites(o.sites.id==siteId,:); % Read the contact axial position.
f = o.patient.fascicles(fascicleId,:); % Choose one provisional fascicle centreline.
%% 2. Interpolate per ampere in the original tetrahedra
z = linspace(-o.params.domain.nerveLength_m/2,o.params.domain.nerveLength_m/2,4001)'; % Sample across the full domain.
v = sampleAxialFields(o.mesh.nodes_m,double(o.mesh.tets),o.V_per_A(:,j),f.x_m,f.y_m,z); % Interpolate in volts per reference current.
assert(all(isfinite(v(2:end-1))),'Missing field inside fibre.'); % Reject internal gaps before coupling.
writematrix([1e6*(z(2:end-1)-site.z_m),v(2:end-1)/o.current_A],outputFile); % Write axial micrometres relative to contact and volts per ampere.
end % Finish exporting one stimulation profile.
