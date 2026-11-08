function af = activatingFunction(potentialFcn, fibres)
%ACTIVATINGFUNCTION Activating function at the nodes of each fibre, for a 1 A stimulus.
%   af = activatingFunction(potentialFcn, fibres)
%
%   Owner: WP4. Status: SKELETON.
%
%   f_n = (V_{n-1} - 2 V_n + V_{n+1}) / dx^2   at each node n of each fibre
%   (dx = internode). f > 0 = depolarisation: a cathodic (negative) AS gives a
%   potential minimum below it -> positive second difference -> excitation.
%
%   Inputs
%     potentialFcn  function handle V = potentialFcn(x_m, y_m, z_m) [V] for 1 A.
%                   From the FEM: @(x,y,z) interpolateSolution(sol.result, x, y, z).
%                   For the sprint-1 test: the analytical point source
%                   @(x,y,z) -1 ./ (4*pi*sigma*sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2)).
%                   (A handle keeps this module independent of the FEM code.)
%     fibres        table from makeFibres
%
%   Output: struct
%     f_Vpm2perA      cell (nFibres,1): f at every node [V/m^2 per A]
%     afMax_Vpm2perA  (nFibres,1) maximum of f along each fibre
%     nodeOfMax       (nFibres,1) index of the node where the maximum occurs

arguments
    potentialFcn (1,1) function_handle
    fibres table
end
error('project:notImplemented', 'activatingFunction: not implemented yet (WP4).');
end
