function values = sampleAxialFields(X,T,U,x,y,z) % Sample any number of nodal fields along one longitudinal fibre.
mask=false(2,2,numel(z)); mask(1,1,:)=true; % Select a single fibre in a narrow regular lattice.
[ids,b]=locateRegularGrid(X,T,{[x,x+1e-6],[y,y+1e-6],z(:)'},mask); % Locate the axial samples inside the true tetrahedra.
ok=isfinite(ids); nodes=T(:,ids(ok))'; % Gather vertices of the containing elements.
values=nan(numel(z),size(U,2)); % Preserve points outside the physical domain.
for k=1:size(U,2) % Interpolate each electrical configuration with the same element map.
    u=U(:,k); values(ok,k)=sum(b(ok,:).*double(u(nodes)),2); % Apply P1 interpolation to the original nodal solution.
end % Finish filling the longitudinal field profiles.
end % Return one row per axial position.
