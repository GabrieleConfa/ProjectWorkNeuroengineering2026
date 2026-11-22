function fp = contactFootprints(sites)
%CONTACTFOOTPRINTS Outline of every distinct contact in the x-y plane (sections +k and -k share it).
fp = struct('corners_m', {}, 'siteIds', {}, 'kind', {});
key = round([sites.centre_m(:,1:2), sites.phi_deg / 1e3] * 1e7);
[~, iu, ic] = unique(key, 'rows', 'stable');
sgn = [-1 -1; 1 -1; 1 1; -1 1];
for k = 1:numel(iu)
    s = sites(iu(k), :);
    fp(k).corners_m = s.centre_m(1:2) + sgn(:,1) * (s.thickness_m/2) * s.n(1:2) ...
                                      + sgn(:,2) * (s.sideXY_m/2) * s.w(1:2);
    fp(k).siteIds = sites.id(ic == k);
    fp(k).kind = string(s.type);
end
end
