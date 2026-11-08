function setup_paths()
%SETUP_PATHS Add the project source and configuration folders to the MATLAB path.
%   Run once per session from anywhere: setup_paths
%   Also warns if the PDE Toolbox (needed for phases 3-4) is not installed.

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'config'));
addpath(fullfile(root, 'scripts'));

if isempty(ver('pde'))
    warning('setup_paths:noPDE', ...
        'Partial Differential Equation Toolbox not found: FEM phases will not run.');
end
end
