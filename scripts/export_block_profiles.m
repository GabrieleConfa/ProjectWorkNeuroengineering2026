%% EXPORT_BLOCK_PROFILES  FEM potential along the fibres, for the conduction-block study (MRG in NEURON).
%
% What it does
%   Reads the solved 3D model results/main3d/<id>_3d.mat (the one shown in the viewer; nothing is
%   re-solved), samples the potential of every solved contact, and of its mirror partner on the other
%   site, along straight fibres inside the fascicles (src/block/sampleFibreProfiles.m), and writes
%     results/block/<id>_fibre_profiles.h5     the profiles (read by scripts/run_block_study.py)
%     results/block/<id>_block_params.json     a copy of config/block_params.m for the Python side
%     results/block/<id>_fibres.png            where the fibres are, over the cross-section
%
% HDF5 layout (MATLAB writes column-major, so h5py sees the dimensions reversed; the shapes below
% are as h5py/numpy read them):
%   /z_m (nZ)   /x_m, /y_m, /fascicle (nL)   /contact_id, /mirrored_from, /contact_z_m (nC)
%   /V_perA (nC, nL, nZ) float32: potential [V] for +1 A injected by that contact alone
%
% Run: setup_paths; export_block_profiles          (then: python scripts/run_block_study.py)
% Time: about a minute (locating ~0.7 M points in the mesh).

clear;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);  setup_paths;
mm = 1e-3;

b = block_params();
S = load(fullfile(root, 'results', 'main3d', b.patientId + "_3d.mat"), 'out');   out = S.out;
outDir = fullfile(root, 'results', 'block');
if ~exist(outDir, 'dir'), mkdir(outDir); end

tic;
prof = sampleFibreProfiles(out, b);
[nL, nZ, nC] = size(prof.V_perA);
fprintf('Fibre profiles: %d fibres x %d samples (|z| <= %.1f mm, every %.0f um) x %d contacts  [%.0f s]\n', ...
        nL, nZ, b.profileHalfSpan_m / mm, b.profileStepZ_m * 1e6, nC, toc);
fprintf('   excluded lines: %d through a contact, %d leaving the endoneurium, %d not found in the mesh\n', ...
        prof.excluded.throughContact, prof.excluded.offEndoneurium, prof.excluded.notFound);
for i = 1:height(out.patient.fascicles)
    fprintf('   F%-2d %3d fibres\n', i, nnz(prof.fascicle == i));
end
sites = out.sites;
cz = arrayfun(@(id) sites.z_m(sites.id == id), prof.contactId);
for j = 1:nC
    src = "solved";
    if prof.mirroredFrom(j) > 0, src = sprintf("mirror of AS %d", prof.mirroredFrom(j)); end
    v = prof.V_perA(:, :, j);
    fprintf('   AS %2d (z = %+5.1f mm, %s): V along the fibres %.1f ... %.1f V/A\n', ...
            prof.contactId(j), cz(j) / mm, src, min(v(:)), max(v(:)));
end

%% Write -------------------------------------------------------------------------------------------
h5 = fullfile(outDir, b.patientId + "_fibre_profiles.h5");
if isfile(h5), delete(h5); end
put = @(name, val) writeH5(h5, name, val);
put('/z_m', prof.z_m(:));
put('/x_m', prof.x_m(:));   put('/y_m', prof.y_m(:));   put('/fascicle', double(prof.fascicle(:)));
put('/contact_id', double(prof.contactId(:)));   put('/mirrored_from', double(prof.mirroredFrom(:)));
put('/contact_z_m', cz(:));
% stored as (nZ, nL, nC) in MATLAB order -> numpy (nC, nL, nZ)
h5create(h5, '/V_perA', [nZ nL nC], 'Datatype', 'single', 'ChunkSize', [nZ 1 1], 'Deflate', 4);
h5write(h5, '/V_perA', permute(prof.V_perA, [2 1 3]));
h5writeatt(h5, '/', 'source', char(prof.source));
h5writeatt(h5, '/', 'created', char(string(datetime("now"))));

js = fullfile(outDir, b.patientId + "_block_params.json");
fid = fopen(js, 'w');   fprintf(fid, '%s', jsonencode(b, PrettyPrint=true));   fclose(fid);
fprintf('Wrote %s\n      %s\n', h5, js);

%% Figure: fibres over the cross-section -------------------------------------------------------------
F = out.patient.fascicles;   ax = out.patient.nerveSemiAxes_m;   rotN = out.patient.nerveRotation_deg;
fig = figure('Visible', 'off', 'Position', [100 100 700 560]);  hold on;  axis equal;
t = linspace(0, 2*pi, 400);
xe = ax(1) * cos(t);   ye = ax(2) * sin(t);
plot((xe * cosd(rotN) - ye * sind(rotN)) / mm, (xe * sind(rotN) + ye * cosd(rotN)) / mm, 'k-');
for i = 1:height(F)
    plot((F.x_m(i) + F.rEndo_m(i) * cos(t)) / mm, (F.y_m(i) + F.rEndo_m(i) * sin(t)) / mm, 'b-');
    text(F.x_m(i) / mm, F.y_m(i) / mm, sprintf('F%d', i), 'Color', 'b', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
end
scatter(prof.x_m / mm, prof.y_m / mm, 6, [0.4 0.4 0.4], 'filled');
st = sites(ismember(sites.id, prof.contactId) & sites.z_m < 0, :);
plot(st.r_m .* cosd(st.theta_deg) / mm, st.r_m .* sind(st.theta_deg) / mm, 'rs', 'MarkerFaceColor', 'r');
text(st.r_m .* cosd(st.theta_deg) / mm, st.r_m .* sind(st.theta_deg) / mm + 0.12, ...
     compose('AS %d', st.id), 'Color', 'r', 'HorizontalAlignment', 'center');
xlabel('x [mm]');  ylabel('y [mm]');
title(sprintf('%s: %d fibre lines (grey), contacts of station z = %.1f mm (red)', b.patientId, nL, ...
      st.z_m(1) / mm));
exportgraphics(fig, fullfile(outDir, b.patientId + "_fibres.png"), 'Resolution', 130);
close(fig);


function writeH5(file, name, val)
%WRITEH5 Create a dataset of the size and class of val and write val into it.
h5create(file, name, size(val), 'Datatype', class(val));
h5write(file, name, val);
end
