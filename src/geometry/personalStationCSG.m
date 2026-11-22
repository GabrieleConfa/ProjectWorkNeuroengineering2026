function [dl,bt]=personalStationCSG(g,Rmed,holes) % Build a contacted cross-section with the original solver's face ordering.
xy=g.outline; n=size(xy,1); f=g.f; nf=size(f,1); S=g.sites; % Read the actual tissue boundary and pad geometry.
gd=zeros(2+2*n,2+2*nf+height(S)); names=strings(1,size(gd,2)); % Allocate MATLAB's CSG primitives.
gd(1:4,1)=[1;0;0;Rmed]; names(1)="M"; % Preserve the grounded outer medium from the 2D benchmark.
gd(:,2)=[2;n;xy(:,1);xy(:,2)]; names(2)="N"; % Define the nerve with the two local flat footprints.
for k=1:nf % Draw each perineurium and endoneurium with their configured semiaxes.
    for inner=0:1 % Preserve the boolean table order expected by both solvers.
        col=2*k+1+inner; axes=f(k,3:4)+(1-inner)*g.q.geometry.tPeri_m; % Add sheath thickness only to the outer ellipse.
        gd(1:6,col)=[4;f(k,1:2)';axes';deg2rad(f(k,5))]; % Encode a rotated ellipse as an exact conic.
        names(col)="F"+k+"L"+inner; % Assign a unique constructive geometry name.
    end % Finish the two tissue interfaces.
end % Finish the fascicles.
signs=[-1 -1;1 -1;1 1;-1 1]; % Use counterclockwise corners for each finite rectangular contact.
for k=1:height(S) % Insert the same oriented volumes visible in the personal viewer.
    col=2+2*nf+k; corners=S.centre_m(k,1:2)+signs(:,1)*S.thickness_m(k)/2*S.n(k,1:2)+signs(:,2)*S.sideXY_m(k)/2*S.w(k,1:2); % Project the box onto this section.
    gd(1:10,col)=[2;4;corners(:,1);corners(:,2)]; names(col)="S"+k; % Keep source contact order.
end % Finish the electrode primitives.
formula=strjoin(names,"+"); % Retain separate metallic regions in the metal workflow.
if holes, formula="("+strjoin(names(1:2+2*nf),"+")+")-("+strjoin(names(3+2*nf:end),"+")+")"; end % Represent ideal electrodes with boundary conditions.
[dl,bt]=decsg(gd,char(formula),char(names(:))'); % Generate the conduction domain.
end
