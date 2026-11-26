function f=midArmLayout(p) % Generate a reproducible synthetic section from histological area targets.
g=p.geometry; spec=g.pw6; n=spec.nFascicles; % Read the shared anatomy targets.
if isfield(g,'layoutTemplate') && size(g.layoutTemplate,1)==n % Use a repeatable reference anatomy for device comparisons.
    template=g.layoutTemplate; centres=template(:,1:2).*g.nerveSemiAxes_m; % Scale the normalized coordinates with both anatomical axes.
    areas=[spec.areaLargest_m2;(spec.areaTotal_m2-spec.areaLargest_m2)*template(2:end,3)/sum(template(2:end,3))]; % Match the two histological area constraints.
    r=sqrt(areas/pi); ang=g.nerveRotation_deg; centres=centres*[cosd(ang) sind(ang);-sind(ang) cosd(ang)]; % Rotate the section as a unit.
    f=[centres r r zeros(n,1)]; return % Use circular equivalent fascicles until individual axes arrive.
end % Other counts use the seeded packing generator.
stream=RandStream('mt19937ar','Seed',g.layoutSeed); % Preserve the user's global random state.
weights=sort(0.5+rand(stream,n-1,1),'descend'); % Allow unequal provisional areas for the remaining fascicles.
areas=[spec.areaLargest_m2;(spec.areaTotal_m2-spec.areaLargest_m2)*weights/sum(weights)]; % Match both the largest and the total area.
r=sqrt(areas/pi); centres=zeros(n,2); % Keep circular fascicles until individual axes are known.
axes=g.nerveSemiAxes_m; angle=(0:127)'*2*pi/128; normals=[cos(angle) sin(angle)]; % Sample the enclosing oval.
support=hypot(axes(1)*normals(:,1),axes(2)*normals(:,2)); % Evaluate the tangent planes of the ellipse.
for attempt=1:250 % Restart if a late fascicle cannot fit.
    for k=1:n % Place larger fascicles first.
        candidates=(2*rand(stream,700,2)-1).*axes; % Sample the section uniformly before excluding overlaps.
        ok=all(candidates*normals'+r(k)+g.tPeri_m+0.035e-3<support',2); % Leave clearance from the epineurial surface.
        for j=1:k-1 % Enforce separation from every previous sheath.
            ok=ok & vecnorm(candidates-centres(j,:),2,2)>r(k)+r(j)+2*g.tPeri_m+0.012e-3; % Preserve a small interfascicular gap.
        end % Finish testing the candidate centres.
        first=find(ok,1); if isempty(first),break;end % Restart this layout when packing fails.
        centres(k,:)=candidates(first,:); % Keep the first admissible random centre.
    end % Finish the current packing.
    if k==n && ~isempty(first) % Return a complete section.
        a=g.nerveRotation_deg; centres=centres*[cosd(a) sind(a);-sind(a) cosd(a)]; % Rotate all centres with the nerve.
        f=[centres r r zeros(n,1)]; return % Return [x y a b angle] in SI.
    end % Try another arrangement.
end % Report incompatible inputs explicitly.
error('Cannot pack the selected fascicle areas within the nerve.');
end
