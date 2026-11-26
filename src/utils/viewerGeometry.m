function model = viewerGeometry(out,xs,ys,zs,gridStepXY_m,gridStepZ_m) % Build the shared geometry metadata with or without electrical fields.
mm=1e-3; p=out.params; patient=out.patient; patientId=patient.id; % Read the actual patient and domain.
sites=out.sites; F=patient.fascicles; R=patient.nerveRadius_m; res=out.results; % Read all inputs required by the viewer.
toMm = @(v) round(v / mm, 5);
solved = out.solvedIds(:)'; % Keep contact ids aligned with the potential columns.
siteList = cell(height(sites), 1);
for k = 1:height(sites)
    j = find(solved == sites.id(k), 1);
    s = struct('id', sites.id(k), 'type', string(sites.type(k)), 'section', sites.section(k), ...
               'centre', toMm(sites.centre_m(k, :)), 'n', round(sites.n(k, :), 6), ...
               'w', round(sites.w(k, :), 6), 'thickness', toMm(sites.thickness_m(k)), ...
               'sideXY', toMm(sites.sideXY_m(k)), 'sideZ', toMm(sites.sideZ_m(k)), ...
               'solved', ~isempty(j), 'location', "", 'R_kOhm', NaN);
    if ~isempty(j)
        s.location = res.location(j);   s.R_kOhm = round(res.R_kOhm(j), 3);
    end
    siteList{k} = s;
end
fascList = cell(height(F), 1);
for i = 1:height(F)
    fascList{i} = struct('id', i, 'x', toMm(F.x_m(i)), 'y', toMm(F.y_m(i)), 'r', toMm(F.rEndo_m(i)), ...
                         'tPeri_um', round(F.tPeri_m(i) * 1e6, 2));
end
model = struct( ...
    'patient', struct('id', patient.id, 'seed', patient.seed, 'level', string(patient.level), ...
                      'nInner', patient.generatorInfo.nInner), ...
    'units', "mm; fields for a 1 uA monopolar stimulus (return = grounded outer boundary)", ...
    'nerveRadius', toMm(R), 'nerveLength', toMm(p.domain.nerveLength_m), ...
    'cuff', struct('length', toMm(p.cuff.length_m), 'wall', toMm(p.cuff.wallThickness_m)), ...
    'outerRadius', toMm(p.domain.outerRadiusFactor * R), ...
    'fascicles', {fascList}, 'sites', {siteList}, 'solvedIds', solved, ...
    'grid', struct('nx', numel(xs), 'ny', numel(ys), 'nz', numel(zs), 'x0', toMm(xs(1)), ...
                   'y0', toMm(ys(1)), 'z0', toMm(zs(1)), 'dx', toMm(gridStepXY_m), ...
                   'dy', toMm(gridStepXY_m), 'dz', toMm(gridStepZ_m)), ...
    'source', sprintf('results/main3d/%s_3d.mat (created %s)', patientId, out.created), ...
    'exported', string(datetime("now")));
if isfield(p,'personalGeometry') % Export the personal schematic overlays from the same input as the contact mesh.
    q = p.personalGeometry; % Read the provisional contact model used for this solve.
    raw = loadActiveSites(q.electrode.sitesFile,R); % Recover the undeviated insertion axes.
    cuff = raw(raw.type=="cuff",:); % Select only the rigid external pads.
    shell = buildQuicklookNerve(q,cuff); % Recalculate the locally indented preview surface.
    ringRows = unique([1;size(shell.x_m,1);find(any(abs(diff(shell.x_m))>1e-12,2));find(any(abs(diff(shell.x_m))>1e-12,2))+1]); % Keep every changing ring and both full-length endpoints.
    for i=1:height(cuff) % Export the contact assumption for each rigid external pad.
        [~,depth,width,measured]=cuffContactPlane(q,cuff(i,:)); % Read the centre displacement and actual contacting width.
        row=find(sites.id==cuff.id(i),1); % Match the physical site instead of its array position.
        model.sites{row}.contact=struct('depth_um',depth*1e6,'area_mm2',width*cuff.sideZ_m(i)*1e6,'measured',measured); % Report finite footprint and provenance in viewer units.
    end % Finish the external pad metadata.
    if isfield(out,'cuffAudit') % Show independently verified mesh triangles only for an audited volume solve.
        elems=double(out.mesh.elemCell); % Map each tetrahedron to its physical region.
        [audit,faces]=checkCuffMeshContact(out.mesh.nodes_m,double(out.mesh.tets),out.cells.tissue(elems),out.cells.site(elems),sites,q); % Reconstruct the actual contact interface.
        for i=1:height(audit) % Attach a small local triangle mesh to each corresponding pad.
            row=find(sites.id==audit.site(i),1); [ids,~,idx]=unique(faces{i}(:)); % Compress the global node indices for this face.
            model.sites{row}.meshContact=struct('area_mm2',audit.shared_mm2(i),'xyz',toMm(out.mesh.nodes_m(:,ids)'),'tri',reshape(idx,size(faces{i}))-1); % Export world coordinates and zero-based indices.
        end % Complete the local interfaces.
    end % Geometry-only previews retain the geometric readout.
    model.nerveAxes = toMm(q.geometry.nerveSemiAxes_m); % Supply the cross-section semiaxes to the viewer.
    model.nerveAngle = q.geometry.nerveRotation_deg; % Supply the ellipse rotation in degrees.
    for k=1:height(F) % Transfer the actual fascicle semiaxes to every rendering mode.
        model.fascicles{k}.axes = toMm(q.geometry.fascicles_m(k,3:4)); % Store width and height in millimetres.
        model.fascicles{k}.angle = q.geometry.fascicles_m(k,5); % Retain the anatomical rotation.
    end % Finish the viewer anatomy.
    model.shellZ = toMm(shell.z_m(ringRows,1)); % Define the ring positions in millimetres.
    model.shellTheta = shell.theta_rad; % Use the actual azimuth of every vertex, including the pad edges.
    model.contactTransition = toMm([q.electrode.cuffLateralTransition_m,q.electrode.cuffTransition_m]); % Supply the same lateral and axial shoulder widths as the MATLAB contact model.
    model.shellR = toMm(hypot(shell.x_m(ringRows,:),shell.y_m(ringRows,:))); % Sample each ring including its repeated closing point.
    model.leadWidth = toMm(q.electrode.stripWidth_m); % Use the width included in the geometric clearance check.
    model.leadThickness = toMm(q.electrode.stripThickness_m); % Display the provisional substrate thickness.
    if isfield(q.geometry,'scenario') && ~isempty(q.geometry.scenario) % Insertion mode: irregular fascicles, PDMS leads and oval PDMS sleeve.
        sc = q.geometry.scenario; model.leadWidth = toMm(sc.leadWidth_m); model.leadThickness = toMm(sc.leadThickness_m); % PDMS lead size.
        zl = unique([-p.domain.nerveLength_m/2; p.domain.nerveLength_m/2; reshape([sc.stations.z_m]+sc.bump_m*[-1;0;1],[],1)]); model.fascZ = toMm(zl(:)'); % Axial levels of the fascicle outlines.
        for k=1:numel(sc.polys0), model.fascicles{k}.poly = cell(1,numel(zl)); end % One outline per level.
        for i=1:numel(zl) % Outlines at every level (displaced near the stations).
            P = scenarioFascicles(sc,zl(i)); for k=1:numel(P), model.fascicles{k}.poly{i} = toMm(P{k}(1:2:end,:)); end % 64 vertices per outline.
        end
        for k=1:numel(sc.polys0), model.fascicles{k}.x = toMm(mean(sc.polys0{k}(:,1))); model.fascicles{k}.y = toMm(mean(sc.polys0{k}(:,2))); end % Label positions: centroids of the undisturbed outlines.
        if isfield(q.geometry,'patient') && ~isempty(q.geometry.patient) % Functions and motor patches of the patient.
            pt = q.geometry.patient; cat = functionCatalog(); model.patientInfo = struct('id',string(pt.id),'source',string(pt.source),'functionSource',string(pt.functionSource)); % Provenance (kept apart: the viewer relies on model.patient.id, which must stay as exported by the solver).
            model.functionCatalog = struct('id',cellstr(cat.id)','name',cellstr(cat.name)','colour',num2cell(cat.colour,2)'); % Group keys, names and colours.
            for k=1:numel(pt.fascicles), c = pt.fascicles(k).composition; model.fascicles{k}.function = ''; model.fascicles{k}.composition = []; % Dominant function (empty = unknown).
                if ~isempty(c), [~,i] = max(c.*(cat.kind=="motor")'); model.fascicles{k}.function = char(cat.id(i)); model.fascicles{k}.composition = round(c,3); end % Dominant motor group (sensory fibres are 87 percent of the axons and would hide it).
                mfk = pt.fascicles(k).motorFraction; if isempty(mfk) || isnan(mfk), mfk = []; else, mfk = round(mfk,3); end, model.fascicles{k}.motorFraction = mfk; model.fascicles{k}.patches = cellfun(@(Q)toMm(Q),pt.fascicles(k).patches,'UniformOutput',false); end
        end
        th = linspace(0,2*pi,97)'; th(end)=[]; ray = [cos(th),sin(th)]; r0 = nerveBoundaryRadius(q,ray); pts = r0(:).*ray; t = circshift(pts,-1,1)-circshift(pts,1,1); nrm = [t(:,2),-t(:,1)]./vecnorm(t,2,2); % Oval sleeve outline.
        model.sleeve = struct('thickness',toMm(sc.sleeveThickness_m),'outer',toMm(pts+sc.sleeveThickness_m*nrm),'inner',toMm(pts),'halfLength',toMm(p.cuff.length_m/2)); % PDMS sleeve around the nerve.
    end
    model.insertionSegments = {}; % Export one continuous centreline at each station.
    for sec = unique(raw.section)' % Process every axial section separately.
        strip = planQuicklookStrip(q,raw(raw.section==sec & raw.type=="intraneural",:),cuff(cuff.section==sec,:),"deflect"); % Recreate the routed contact path.
        model.insertionSegments{end+1}=toMm(strip.path_m); % One connected support per station, terminating after its fourth pad.
    end % Finish all eight sections.
end % Reference patients retain their original cylindrical shell.
end % Return the geometry for export.
