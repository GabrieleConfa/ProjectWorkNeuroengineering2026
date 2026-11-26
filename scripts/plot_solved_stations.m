function plot_solved_stations(matFile) % Draw six solved contact fields on each actual anatomical section.
S=load(matFile,'out'); o=S.out; % Read the saved mesh and all contact potentials.
q=o.params.personalGeometry; % Recover geometry used by the full solver.
raw=loadActiveSites(q.electrode.sitesFile,q.geometry.radius_m); % Read the four nominal pads and two exterior contacts per station.
xx=linspace(-o.patient.nerveRadius_m,o.patient.nerveRadius_m,401); % Sample the transverse field in metres (about 11 µm per pixel, so the map and the nerve edge look smooth).
[GX,GY]=ndgrid(xx,xx); % Use x as the fastest varying display index.
angles=linspace(0,2*pi,361)'; ray=[cos(angles),sin(angles)]; % Sample the actual nerve circumference.
dirName=fileparts(matFile); % Write the plots alongside the numerical solution.
for sec=unique(raw.section)' % Draw all eight axial levels independently.
    rows=find(raw.section==sec); z=raw.z_m(rows(1)); % Locate this exact plane from the tutor's coordinates.
    cuff=raw(rows(raw.type(rows)=="cuff"),:); % Select the two rigid exterior contacts.
    strip=planQuicklookStrip(q,raw(rows(raw.type(rows)=="intraneural"),:),cuff,"deflect"); % Recover the physical lead.
    rad=cuffDeformedRadius(q,cuff,z,ray); outline=rad(:).*ray; % Calculate local apposition on this section.
    mask=inpolygon(GX,GY,outline(:,1),outline(:,2)); % Retain only nerve tissue.
    both=cat(3,mask,false(size(mask))); % Use a two-plane sampling lattice with only the actual station selected.
    [ids,b]=locateRegularGrid(o.mesh.nodes_m,double(o.mesh.tets),{xx,xx,[z,z+1e-3]},both); % Locate each transverse sample in the full mesh.
    hit=find(mask); ok=isfinite(ids); hit=hit(ok); % Exclude points on the faceted exterior missed by the mesh.
    nodes=double(o.mesh.tets(:,ids(ok))'); % Gather the four interpolation nodes for each sample.
    fig=figure(Visible="off",Theme="light",Color="white",Position=[60 60 1250 760]); % Create a readable six-contact sheet.
    layout=tiledlayout(fig,2,3,TileSpacing="compact"); % Show all six stimulus choices.
    for j=1:numel(rows) % Plot one monopolar solution per panel.
        i=rows(j); col=find(o.solvedIds==raw.id(i),1); % Look up the correct solution by contact ID.
        u=double(o.V_per_A(:,col))*1e-3/o.current_A; % Scale the reference to mV for 1 µA.
        v=nan(size(mask)); v(hit)=sum(b(ok,:).*u(nodes),2); % Interpolate voltage in actual tissue tetrahedra.
        ax=nexttile(layout); h=imagesc(ax,xx*1e3,xx*1e3,v'); % Display physical coordinates in millimetres.
        h.AlphaData=isfinite(v'); axis(ax,'xy','equal','tight'); hold(ax,'on'); % Hide exterior and preserve proportions.
        cb=colorbar(ax); cb.Label.String='mV'; % Label the per-panel potential scale.
        plot(ax,outline(:,1)*1e3,outline(:,2)*1e3,'Color',[0.1 0.3 0.85]); % Show the finite cuff fit.
        for P=stationFascicles(q,sec)', plot(ax,[P{1}(:,1);P{1}(1,1)]*1e3,[P{1}(:,2);P{1}(1,2)]*1e3,'Color',[0.45 0.45 0.45]); end % Fascicle outlines of this station (displaced by its lead).
        plot(ax,strip.path_m(:,1)*1e3,strip.path_m(:,2)*1e3,'Color',[0.12 0.36 0.94],'LineWidth',1.7); % Draw the single-entry carrier.
        for k=rows' % Draw the six active faces on their local supports.
            s=o.sites(k,:); centre=s.centre_m(1:2); % Read the FEM contact centre.
            if s.type=="cuff", centre=centre-s.thickness_m/2*s.n(1:2); end % Use the inner active cuff face.
            edge=centre+[-1;1]*s.sideXY_m/2*s.w(1:2); % Compute transverse active-face endpoints.
            plot(ax,edge(:,1)*1e3,edge(:,2)*1e3,'Color',[0.62 0.13 0.85],'LineWidth',3); % Display all pads in purple.
        end % Complete the device overlay.
        title(ax,sprintf('AS %d | %s | %.2f kΩ',raw.id(i),string(raw.type(i)),o.results.R_kOhm(col))); % Identify the energized contact.
        xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); % Label both spatial axes.
    end % End the six electrical configurations.
    title(layout,sprintf('Station %+d | z = %.1f mm | monopolar 1 µA',sec,z*1e3)); % Identify the physical section.
    exportgraphics(fig,fullfile(dirName,sprintf('Cuff60_field_%+d.png',sec)),Resolution=200); close(fig); % Save the final section and release its graphics.
end % Finish the complete station atlas.
end % Return after saving every field sheet.
