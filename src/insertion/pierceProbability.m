function P = pierceProbability(r_m, gamma_deg, ip)
%PIERCEPROBABILITY Probability that the lead pierces a fascicle it touches (pierceModel = "probability").
%   P = pierceProbability(r_m, gamma_deg, ip)
%
%   Purpose
%     The tutor (meeting 2026-11-23) described the insertion as: at each fascicle on its way the probe
%     pierces it with a probability that grows with the fascicle size, the small fascicles of our
%     model having a very low probability and the large ones a high one; otherwise the probe bends
%     around it. The user (2026-11-26) asked for the probability to depend also on the angle of
%     penetration. The law is the product of two factors:
%
%       P = P_r(r) * P_gamma(gamma)
%       P_r(r)         = 1 / (1 + (r50 / r)^n)     size factor: a smooth step from 0 (small) to 1 (large);
%                                                  r50 = radius with P_r = 50 %, n = steepness
%       P_gamma(gamma) = cos(gamma)^m              angle factor: 1 head-on (gamma = 0), 0 tangent (90 deg);
%                                                  only the normal component of the push, F cos(gamma),
%                                                  loads the perineurium (the same F cos(gamma) as in
%                                                  Marco's force model); m sets how fast a grazing
%                                                  contact turns into a deflection
%
%     With the default values (config/insertionParams.m: r50 = 0.38 mm, n = 6, m = 2), on patient UA01
%     (fascicle radii 0.17 ... 0.69 mm) a head-on contact pierces F10 (0.17 mm) with ~1 %, F7
%     (0.28 mm) ~13 %, F5 (0.34 mm) ~34 %, F3 (0.43 mm) ~68 %, F2 (0.57 mm) ~92 %, F1 (0.69 mm)
%     ~97 %; at 45 deg these are halved, at 60 deg divided by 4. All three values are [ASSUMED],
%     chosen to reproduce the tutor's qualitative description: the law and its values are to be asked
%     to the tutor (or fitted to insertion data, e.g. the angle dependence reported by Sergi, Jensen &
%     Yoshida, arXiv:1308.0525, not yet read).
%
%   Inputs
%     r_m        equivalent radius sqrt(A/pi) of the fascicle [m] (scalar or array)
%     gamma_deg  angle between the lead heading and the outward normal of the perineurium at the
%                contact [deg]; 0 = head-on (same size as r_m, or scalar)
%     ip         insertion parameters (pierceR50_m, pierceSteepness, pierceAngleExponent)
%
%   Output
%     P          probability in [0, 1]

arguments
    r_m double {mustBeNonnegative}
    gamma_deg double
    ip (1,1) struct
end
Pr = 1 ./ (1 + (ip.pierceR50_m ./ max(r_m, eps)).^ip.pierceSteepness);
Pg = max(0, cosd(gamma_deg)).^ip.pierceAngleExponent;
P = Pr .* Pg;
end
