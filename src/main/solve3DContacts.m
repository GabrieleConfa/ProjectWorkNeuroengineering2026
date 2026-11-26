function [V, res, key] = solve3DContacts(M, solveIds, opts)
%SOLVE3DCONTACTS Solve the 3D model for a unit current on each of the given contacts.
%   [V, res, key] = solve3DContacts(M, solveIds)
%   [V, res, key] = solve3DContacts(M, solveIds, opts)
%
%   Purpose
%     Section 4 (solve) of main_pipeline_3d as a function (split out on 2026-11-26 without changing
%     the computation). For each contact: the contact becomes the conductor (idle contacts
%     insulating, D-07), a uniform source I0/volume is spread over its elements, and the system
%     K u = F is solved with the grounded faces at 0 V. Linear problem: the solution for any current
%     I is I/I0 times it.
%
%   Per-contact cache (opts.cacheDir). Every contact is written to disk as soon as it is solved
%   (cacheDir/AS<nn>.mat), so that a long run that stops (out of memory, MATLAB closed) can be
%   resumed: the contacts already on disk are loaded instead of solved. A file is reused only if it
%   was computed on the SAME FE problem: same mesh and same solve parameters (solveKey below). The
%   mesh is regenerated on every run, so resuming relies on generateMesh giving the same mesh for the
%   same geometry (it does: checked on 2026-11-26); if not, the contacts are solved again.
%
%   Inputs
%     M         struct from build3DModel
%     solveIds  contact ids to solve
%     opts      optional struct:
%       .cacheDir  ("")  folder of the per-contact cache; "" = no cache (studies)
%       .prevFile  ("")  a complete earlier result (results/main3d/<id>_3d.mat) whose contacts seed the
%                        cache if it was computed on the same FE problem
%       .verbose   (true)
%
%   Outputs
%     V    nAll x nS potentials [V per I0] (single if p.solve3d.saveSingle)
%     res  table per contact: id, location, Vcontact_V, R_kOhm, spread_pct, conservErr_pct
%     key  the fingerprint of the FE problem (solveKey)

arguments
    M (1,1) struct
    solveIds (:,1) double
    opts (1,1) struct = struct()
end
if ~isfield(opts, 'cacheDir'), opts.cacheDir = ""; end
if ~isfield(opts, 'prevFile'), opts.prevFile = ""; end
if ~isfield(opts, 'verbose'),  opts.verbose = true; end
mm = 1e-3;
p = M.p;   sites = M.sites;   tl = M.tl;   nAll = M.nAll;
nS = numel(solveIds);   I0 = p.stim.unitCurrent_A;
if p.solve3d.saveSingle, V = zeros(nAll, nS, 'single'); else, V = zeros(nAll, nS); end
res = table(solveIds, strings(nS,1), zeros(nS,1), zeros(nS,1), zeros(nS,1), zeros(nS,1), ...
            'VariableNames', {'id', 'location', 'Vcontact_V', 'R_kOhm', 'spread_pct', 'conservErr_pct'});
key = solveKey(tl.nodes_m, tl.tets, p);
useCache = strlength(string(opts.cacheDir)) > 0;
if useCache && ~exist(opts.cacheDir, 'dir'), mkdir(opts.cacheDir); end

% Seed the cache from a previous complete run, if it was computed on the same FE problem.
if useCache && strlength(string(opts.prevFile)) > 0 && isfile(opts.prevFile)
    prev = load(opts.prevFile, 'out');   prev = prev.out;
    if isequaln(solveKey(prev.mesh.nodes_m, prev.mesh.tets, prev.params), key)
        for jj = 1:numel(prev.solvedIds)
            cfile = fullfile(opts.cacheDir, sprintf('AS%02d.mat', prev.solvedIds(jj)));
            if ~isfile(cfile)
                u_save = prev.V_per_A(:, jj);   resRow = prev.results(jj, :);   %#ok<NASGU>
                save(cfile, 'u_save', 'resRow', 'key', '-v7');
            end
        end
        fprintf('\nCache seeded with the %d contacts of %s (same mesh and parameters)\n', ...
                numel(prev.solvedIds), opts.prevFile);
    else
        fprintf('\n%s was computed on a different mesh or parameters: its contacts are not reused\n', opts.prevFile);
    end
    clear prev
end

if opts.verbose
    fprintf('\nSolving %d contacts (unit current %g A, idle contacts %s):\n', nS, I0, p.electrode.idleContacts);
end
for j = 1:nS
    s = solveIds(j);   row = find(sites.id == s);
    cfile = "";
    if useCache
        cfile = fullfile(opts.cacheDir, sprintf('AS%02d.mat', s));
        if isfile(cfile)                                   % already solved by an earlier run?
            C = load(cfile);
            if isequaln(C.key, key)
                V(:, j) = C.u_save;   res(j, :) = C.resRow;
                fprintf('   AS %2d loaded from %s\n', s, cfile);
                continue
            end
            fprintf('   AS %2d: the file on disk is for another mesh or parameters, solving again\n', s);
        end
    end
    es = find(M.elemSite == s);
    K = M.K0;
    if p.electrode.idleContacts == "insulating"           % this contact becomes the conductor
        K = K + assembleConductionP1(tl.nodes_m, tl.tets(:, es), ...
                                     (p.electrode.contactSigma_Spm - M.sIdle) * [1; 1; 1]);
    end
    Te = tl.tets(:, es);
    F = accumarray(Te(:), reshape(repmat(I0 / sum(M.vol(es)) * M.vol(es) / 4, 4, 1), [], 1), [nAll 1]);
    u = zeros(nAll, 1);   tic;
    u(M.free) = solveSPD(K(M.free, M.free), F(M.free), p.solve3d);
    Ignd = -sum(K(M.gnd, :) * u - F(M.gnd));               % current leaving through the ground
    cn = unique(Te(:));   Vc = mean(u(cn));
    res.location(j) = M.where(row);   res.Vcontact_V(j) = Vc;   res.R_kOhm(j) = Vc / I0 / 1e3;
    res.spread_pct(j) = 100 * (max(u(cn)) - min(u(cn))) / Vc;
    res.conservErr_pct(j) = 100 * abs(Ignd - I0) / I0;
    V(:, j) = u;
    if opts.verbose
        fprintf('   AS %2d (%s, z = %+.1f mm) %-28s R = %7.2f kOhm, contact spread %.3f%%, current check %.1e%%  [%.0f s]\n', ...
                s, string(sites.type(row)), sites.z_m(row)/mm, M.where(row), res.R_kOhm(j), ...
                res.spread_pct(j), res.conservErr_pct(j), toc);
    end
    if useCache                                            % write this contact to disk at once
        u_save = V(:, j);   resRow = res(j, :);            %#ok<NASGU>
        save(cfile, 'u_save', 'resRow', 'key', '-v7');
    end
end
end


function key = solveKey(nodes_m, tets, p)
%SOLVEKEY Fingerprint of the FE problem, used to decide whether a contact solved by an earlier run
%   can be reused: it changes if the mesh (number of nodes and tetrahedra, node coordinates,
%   connectivity) or any parameter that enters the solve changes.
key = struct('nNodes', size(nodes_m, 2), 'nTets', size(tets, 2), ...
             'nodeSum', sum(double(nodes_m), 'all'), 'nodeAbsSum', sum(abs(double(nodes_m)), 'all'), ...
             'tetSum', sum(double(tets), 'all'), ...
             'materials', p.materials, 'perineurium', p.perineurium, 'electrode', p.electrode, ...
             'current_A', p.stim.unitCurrent_A, 'solver', p.solve3d.solver);
end


function x = solveSPD(A, b, opt)
%SOLVESPD Solve the symmetric positive definite system A x = b (direct or preconditioned CG).
switch opt.solver
    case "direct"
        x = A \ b;
    case "pcg"
        % Conjugate gradients with an incomplete-Cholesky preconditioner (default since 2026-11-26,
        % D-14): much less memory than the direct factorisation, which nearly filled 32 GB. A
        % non-converged solution must not be saved as a result: a stagnation just above the
        % tolerance is accepted below 1e-8 (Marco's branch); otherwise, fall back to the direct
        % solver. The current balance (conservErr_pct) is checked afterwards anyway.
        Lc = ichol(A, struct('type', 'ict', 'droptol', 1e-4, 'diagcomp', 1e-3));
        [x, flag, relres, iters] = pcg(A, b, opt.pcgTol, 5000, Lc, Lc');
        if ~(flag == 0 || relres < 1e-8)
            warning('solve3DContacts:pcg', 'pcg did not converge (flag %d, relres %.1e after %d iterations): direct solve instead.', ...
                    flag, relres, iters);
            clear Lc
            x = A \ b;
        end
    otherwise
        error('Unknown solver "%s".', opt.solver);
end
end
