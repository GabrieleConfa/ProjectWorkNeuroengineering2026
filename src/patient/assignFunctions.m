function patient = assignFunctions(patient, scenario, seed, opts) % Give every fascicle of a VIRTUAL patient a random functional composition.
%% Purpose
% At mid-arm there is no published map from fascicle to muscle (Delgado-Martinez 2016: "no prediction of destination can be
% established a priori"; Sunderland 1945 found no consistent pattern). The assignment is therefore a stated hypothesis. Three
% hypotheses are available; a measured patient does not use this function (its composition comes from the data file).
%   "random"      (scenario C) every fascicle gets an independent random composition (Dirichlet) around the prior group weights.
%   "forearmRule" (scenario A) the forearm rule extended proximally: forearm muscle groups are more frequent in peripheral and
%                  dorsal fascicles, hand groups in central ones. Needs the orientation of the section: ASSUMED +y = dorsal.
%   "mrAIN"       (scenario B) MR neurography: the anterior interosseous group sits in dorsal and lateral fascicles; the other
%                  groups are random. ASSUMED orientation: +y = dorsal, +x = lateral (radial).
% Motor fraction of every fascicle: about 13 percent on average (Delgado-Martinez 2016: 14 and 12 percent), drawn from a Beta law.
if nargin<4, opts = struct(); end
if ~isfield(opts,'concentration'), opts.concentration = 3; end % ASSUMED: how strongly the groups are mixed inside a fascicle (small = mixed, large = close to the prior).
if ~isfield(opts,'strength'), opts.strength = 0.8; end % ASSUMED: strength of the spatial rules A and B (0 = none, 1 = strong).
if ~isfield(opts,'motorMean'), opts.motorMean = 0.13; end % Mean motor fraction of the axons (literature).
cat = functionCatalog(); rs = RandStream('mt19937ar','Seed',seed); n = numel(patient.fascicles);
isMotor = cat.kind=="motor"; wM = cat.priorWeight(isMotor)/sum(cat.priorWeight(isMotor)); wS = cat.priorWeight(~isMotor)/sum(cat.priorWeight(~isMotor));
C = cell2mat(arrayfun(@(f)mean(f.outline,1),patient.fascicles(:),'UniformOutput',false)); % Fascicle centroids (m).
nerveC = mean(patient.nerve.outline,1); rel = C-nerveC; rmax = max(vecnorm(patient.nerve.outline-nerveC,2,2));
rho = min(1,vecnorm(rel,2,2)/rmax); dorsal = rel(:,2)/rmax; lateral = rel(:,1)/rmax; % Peripheral index (0 centre, 1 border), dorsal and lateral coordinates.
kM = 6; mf = betaRand(rs,opts.motorMean*kM,(1-opts.motorMean)*kM,n); % Motor fraction per fascicle (Beta, mean 0.13).
for k = 1:n
    wm = wM; s = opts.strength; % Group weights of the motor fibres in this fascicle.
    switch scenario
        case "forearmRule"
            forearm = ismember(cat.id(isMotor),["PT","FCRPL","FDS","AIN"]); peripheral = 0.5*rho(k)+0.25*(dorsal(k)+1); % 0 to 1: peripheral and dorsal.
            wm = wm.*(forearm.*((1-s)+s*peripheral)+~forearm.*((1-s)+s*(1-peripheral)));
        case "mrAIN"
            isA = cat.id(isMotor)=="AIN"; score = 0.5*(dorsal(k)+lateral(k)); wm(isA) = wm(isA)*exp(3*s*score); % AIN fibres favour dorsal and lateral fascicles.
        case "random"
        otherwise, error('Unknown scenario "%s".',scenario);
    end
    wm = wm/sum(wm); comp = zeros(height(cat),1);
    comp(isMotor) = mf(k)*dirichletRand(rs,opts.concentration*wm); % Motor part of the composition.
    comp(~isMotor) = (1-mf(k))*dirichletRand(rs,opts.concentration*wS); % Sensory part.
    patient.fascicles(k).composition = comp'; patient.fascicles(k).motorFraction = mf(k); % Store the composition (same order as the catalogue).
    patient.fascicles(k).patches = motorPatches(patient.fascicles(k).outline,mf(k),rs); % Motor patches inside the fascicle.
end
patient.functionScenario = char(scenario); patient.functionSource = 'assumed'; patient.functionIds = cellstr(cat.id)';
end
function x = betaRand(rs,a,b,n) % Beta random numbers from two gamma variables (no toolbox).
x = zeros(n,1); for i = 1:n, g1 = gammaRand(rs,a); g2 = gammaRand(rs,b); x(i) = g1/(g1+g2); end
end
function d = dirichletRand(rs,alpha) % Dirichlet random vector.
g = arrayfun(@(a)gammaRand(rs,a),alpha(:)); d = g/sum(g);
end
function g = gammaRand(rs,a) % Gamma(a,1) by Marsaglia and Tsang (a > 0).
if a<1, g = gammaRand(rs,a+1)*rand(rs)^(1/a); return; end
d = a-1/3; c = 1/sqrt(9*d);
while true
    x = randn(rs); v = (1+c*x)^3; if v<=0, continue; end
    u = rand(rs); if log(u) < 0.5*x^2+d-d*v+d*log(v), g = d*v; return; end
end
end
