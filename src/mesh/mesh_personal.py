"""Conforming linear tetrahedra for the finite rigid contacts. CAD units are mm."""
import json  # Read the MATLAB device description.
import math  # Rotate CAD coordinates in radians.
import os  # Read the optional diagnostic switch.
import sys  # Read input and output paths from MATLAB.
import gmsh  # Use the OpenCASCADE and mesh APIs.
import numpy as np  # Store nodal coordinates and tetrahedral connectivity.
from scipy.io import savemat  # Return the mesh as a MATLAB data file.

def poly_wire(occ, xy, z):  # Closed polygon wire at height z.
    points = [occ.addPoint(x, y, z) for x, y in xy]  # Vertices.
    return occ.addWire([occ.addLine(points[i], points[(i+1)%len(points)]) for i in range(len(points))])  # Closed wire.

def build(cfg):  # Build and label the shared contact mesh.
    # 1. Nerve, fascicles, muscle and rigid pad solids
    occ = gmsh.model.occ  # Select the solid geometry kernel.
    length, cuff, radius, wall = cfg["domain"]  # Read the axial lengths and surrounding radii.
    axial_scale = 0.5  # Stretch tetrahedra along the longitudinal fibre direction after meshing.
    length *= axial_scale  # Scale the full domain for longitudinal mesh spacing.
    cuff *= axial_scale  # Scale the cuff length by the same factor.
    z0 = -length / 2  # Centre the device on the axial origin.
    def prism(face):  # Extrude a fascicle profile across the physical domain.
        return next(tag for dim, tag in occ.extrude([(2, face)], 0, 0, length) if dim == 3)
    levels = np.array(cfg["sectionZ"]) * axial_scale  # Convert station positions into meshing coordinates.
    directions = np.asarray(cfg["rays"])  # Match each vertex with its counterpart on the next section.
    wires = []  # Store the ordered wires of the solved nerve.
    for level, radii in zip(levels, cfg["radii"]):  # Start with the proximal end and finish at the distal end.
        points = [occ.addPoint(x, y, level) for x,y in directions * np.asarray(radii)[:,None]]  # Place each deformed ring.
        edges = [occ.addLine(point, points[(i+1)%len(points)]) for i,point in enumerate(points)]  # Close the transverse profile.
        wires.append(occ.addWire(edges))  # Preserve consistent winding across stations.
    loft = occ.addThruSections(wires, makeSolid=True, makeRuled=True)  # Extrude matching profiles and linearly connect differing sections before meshing.
    bodies = [(next(tag for dim,tag in loft if dim==3), (3,0,0))]  # Assign the connected solid to epineurium.
    if "sleeveXY" in cfg:  # Insertion mode: PDMS sleeve and full muscle cylinder.
        bodies.append((occ.addCylinder(0, 0, z0, 0, 0, length, radius), (4, 0, 0)))  # Muscle around the whole domain (the nerve and sleeve take priority).
        sleeve = [poly_wire(occ, xy, zz*axial_scale) for xy, zz in zip(cfg["sleeveXY"], cfg["sleeveZ"])]  # Offset contours along the cuff length.
        sl = occ.addThruSections(sleeve, makeSolid=True, makeRuled=True)  # Oval PDMS sleeve solid (it contains the nerve; the nerve wins the overlap).
        bodies.append((next(tag for dim, tag in sl if dim == 3), (7, 0, 0)))  # Label as sleeve.
    else:  # Legacy: circular insulating void around the nerve.
        ring = occ.addCylinder(0, 0, -cuff/2, 0, 0, cuff, radius)  # Create the surrounding muscle blank.
        bore = occ.addCylinder(0, 0, -cuff/2, 0, 0, cuff, wall)  # Represent the insulating cuff region.
        annulus, _ = occ.cut([(3, ring)], [(3, bore)])  # Remove that region from the muscle.
        bodies.append((annulus[0][1], (4, 0, 0)))  # Keep the conducting annulus.
        for start in [z0, cuff/2]:  # Build the two muscle end regions.
            bodies.append((occ.addCylinder(0, 0, start, 0, 0, (length-cuff)/2, radius), (4, 0, 0)))  # Label muscle beyond the cuff.
    if "fascXY" in cfg:  # Irregular fascicles: ruled loft through their outlines at every axial level.
        for k, (sections, levels) in enumerate(zip(cfg["fascXY"], cfg["fascZ"]), 1):
            wires = [poly_wire(occ, xy, zz*axial_scale) for xy, zz in zip(sections, levels)]  # Outline at each level.
            fl = occ.addThruSections(wires, makeSolid=True, makeRuled=True)  # Fascicle solid.
            bodies.append((next(tag for dim, tag in fl if dim == 3), (1, k, 0)))  # Endoneurium and fascicle label.
    for k, (x, y, a, b, angle) in enumerate([] if "fascXY" in cfg else cfg["fascicles"], 1):  # Legacy ellipse fascicles.
        # Four exact ellipse arcs keep each curved patch shorter than a semicircle.
        points = []  # Collect the four rotated ellipse vertices.
        for theta in [0, math.pi/2, math.pi, 3*math.pi/2]:  # Visit each end of the two principal axes.
            u, v = a*math.cos(theta), b*math.sin(theta)  # Locate the vertex in the local frame.
            phi = math.radians(angle)  # Convert the supplied orientation to radians.
            points.append(occ.addPoint(x+u*math.cos(phi)-v*math.sin(phi), y+u*math.sin(phi)+v*math.cos(phi), z0))  # Place the vertex in the anatomical frame.
        center = occ.addPoint(x,y,z0)  # Mark the centre used to define ellipse arcs.
        major = points[0 if a >= b else 1]  # Identify the major-axis endpoint for OpenCASCADE.
        arcs = [occ.addEllipseArc(point,center,major,points[(i+1)%4]) for i,point in enumerate(points)]  # Join the four exact quarter arcs.
        bodies.append((prism(occ.addPlaneSurface([occ.addWire(arcs)])), (1, k, 0)))  # Assign the endoneurium and fascicle labels.
    for lead in cfg.get("leads", []):  # PDMS body of each lead: polygon band extruded over its axial thickness.
        w = poly_wire(occ, lead["xy"], (lead["z"]-lead["thickness"]/2)*axial_scale)  # Band outline at the lower face.
        face = occ.addPlaneSurface([w])  # Planar face.
        ext = occ.extrude([(2, face)], 0, 0, lead["thickness"]*axial_scale)  # Extrude to the upper face.
        bodies.append((next(tag for dim, tag in ext if dim == 3), (6, 0, 0)))  # PDMS lead body.
    for site, x, y, z, thick, width, height, phi in cfg["pads"]:  # Build each finite electrode in its own frame.
        z *= axial_scale  # Scale its axial centre coordinate.
        height *= axial_scale  # Scale its axial size consistently.
        box = occ.addBox(-thick/2, -width/2, -height/2, thick, width, height)  # Centre a rigid contact at the local origin.
        occ.rotate([(3, box)], 0, 0, 0, 0, 0, 1, math.radians(phi))  # Align its normal with the MATLAB pose.
        occ.translate([(3, box)], x, y, z)  # Move it to the actual pad centre.
        bodies.append((box, (5, 0, int(site))))  # Retain the electrode identifier.

    # 2. Conforming material interfaces
    inputs = [(3, tag) for tag, _ in bodies]  # List the solids participating in the same mesh.
    _, pieces = occ.fragment(inputs[:1], inputs[1:])  # Split solids at material interfaces.
    occ.synchronize()  # Publish the shared boundaries to the mesher.
    priority = {4: 0, 7: 1, 3: 2, 1: 3, 6: 4, 5: 5}  # Give the inner material ownership of intersecting cells.
    roles = {}  # Map fragmented cells back to physical materials.
    for (_, role), fragments in zip(bodies, pieces):  # Read each original solid and its fragments.
        for dim, tag in fragments:  # Visit each fragment of this solid.
            if dim == 3 and (tag not in roles or priority[role[0]] > priority[roles[tag][0]]):  # Resolve overlapping source labels on volumes.
                roles[tag] = role  # Record the winning tissue and device labels.
    # 3. Electrode and tissue mesh sizes
    field = gmsh.model.mesh.field  # Create spatial mesh-size controls.
    surfaces = set()  # Collect the physical electrode surfaces.
    for tag, role in roles.items():  # Read the labels of the final volumes.
        if role[0] == 5:  # Select the conducting metal.
            surfaces.update(t for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False))  # Include every face of this electrode.
    contact, endo, nerve, coarse, maximum, active = cfg["sizes"]  # Read element sizes and the active half-length.
    distance = field.add("Distance")  # Measure distance from the active metal.
    field.setNumbers(distance, "SurfacesList", sorted(surfaces))  # Use the union of the physical pad faces.
    threshold = field.add("Threshold")  # Transition smoothly to coarser elements.
    field.setNumber(threshold, "InField", distance)  # Use electrode distance to control refinement.
    for name, val in dict(SizeMin=contact, SizeMax=maximum, DistMin=contact, DistMax=3).items():  # Set near and far size limits in millimetres.
        field.setNumber(threshold, name, val)  # Apply that limit to this field.
    radial = field.add("MathEval")  # Add the independent tissue-size control.
    nerve_r = max(math.hypot(x,y) for x,y in cfg["outline"])  # Measure the largest native radius.
    field.setString(radial, "F", f"Min({maximum}, {endo}+{nerve-endo}*(Min(1,Sqrt(x*x+y*y)/{nerve_r}))+0.4*Max(0,Sqrt(x*x+y*y)-{nerve_r})+Min({coarse-nerve},0.3*Max(0,Abs(z)/{axial_scale}-{active})))")  # Grade the element size radially and along the nerve.
    fields = [threshold, radial]  # Size controls that are always active.
    edge_size = cfg.get("edgeSize", 0)  # Optional element size on the pad edges in mm (0 = off).
    if edge_size:  # The current concentrates on pad edges and corners, so resolve them locally.
        wanted = set(cfg.get("edgePads", []))  # Pads whose edges are refined (empty = all pads).
        edge_faces = set()  # Faces of the selected pads.
        for tag, role in roles.items():  # Visit every labelled volume.
            if role[0] == 5 and (not wanted or int(role[2]) in wanted):  # A selected electrode.
                edge_faces.update(t for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False))  # Its faces.
        curves = sorted({t for d, t in gmsh.model.getBoundary([(2, s) for s in edge_faces], combined=False, oriented=False)})  # Every edge of the selected pads.
        edge_distance = field.add("Distance")  # Measure distance from the pad edges.
        field.setNumbers(edge_distance, "CurvesList", curves)  # Use the union of the pad edges.
        field.setNumber(edge_distance, "Sampling", 20)  # Sample each edge densely enough for small elements.
        edge_threshold = field.add("Threshold")  # Grow from the edge size to the coarse size over 0.15 mm.
        field.setNumber(edge_threshold, "InField", edge_distance)
        for name, val in dict(SizeMin=edge_size, SizeMax=maximum, DistMin=edge_size, DistMax=0.15).items():  # SizeMax applies everywhere beyond DistMax, so it must be the coarsest size.
            field.setNumber(edge_threshold, name, val)
        fields.append(edge_threshold)  # Add it to the combined size field.
        if os.environ.get("NERVE_DEBUG_EDGES"):  # Diagnostic: report the selected edge length and stop before meshing.
            total = sum(occ.getMass(1, c) for c in curves); print(f"EDGE CURVES {len(curves)} total length {total:.3f} mm", flush=True); sys.exit(0)
    minimum = field.add("Min")  # Choose the finer local size everywhere.
    field.setNumbers(minimum, "FieldsList", fields)  # Combine electrode, edge and tissue refinement.
    field.setAsBackgroundMesh(minimum)  # Activate the combined size field.
    for option, number in {"General.NumThreads":4,"Mesh.Algorithm3D":1,"Mesh.OptimizeThreshold":0.1,"Mesh.MeshSizeExtendFromBoundary":0,"Mesh.MeshSizeFromPoints":0,"Mesh.MeshSizeFromCurvature":0}.items():  # Set the tetrahedral algorithm and size controls.
        gmsh.option.setNumber(option, number)  # Apply the current meshing option.
    import time; print('CAD done', flush=True)  # Progress marker.
    gmsh.model.mesh.generate(3)  # Tetrahedralize all conforming cells together.

    # 4. Physical coordinates, connectivity and signed-volume checks
    tags, xyz, _ = gmsh.model.mesh.getNodes()  # Get the global coordinates and node identifiers.
    xyz = xyz.reshape(-1,3)  # Store one physical point per row.
    xyz[:,2] /= axial_scale  # Restore the true device length before assembly.
    order = np.argsort(tags)  # Sort the original Gmsh node identifiers.
    sorted_tags = tags[order]  # Use sorted indices to build MATLAB connectivity.
    tets, cells = [], []  # Collect tetrahedra and their parent cells.
    for i, tag in enumerate(roles, 1):  # Visit all labelled volume cells.
        types, _, nodes = gmsh.model.mesh.getElements(3, tag)  # Read the connectivity belonging to this cell.
        assert list(types) == [4], f"Expected linear tetrahedra in volume {tag}"  # Require linear tetrahedra for the P1 solver.
        elements = np.searchsorted(sorted_tags, nodes[0]).reshape(-1, 4)+1  # Convert node identifiers to one-based indices.
        tets.append(elements)  # Keep the block in the merged connectivity.
        cells.append(np.full(len(elements), i, dtype=np.uint32))  # Label each tetrahedron with its cell index.
    all_tets = np.vstack(tets)  # Concatenate the cell blocks.
    coordinates = xyz[order]  # Match coordinates to sorted node identifiers.
    for start in range(0, len(all_tets), 50000):  # Check signed volumes in bounded blocks.
        vertices = coordinates[all_tets[start:start+50000]-1]  # Gather four corners per physical tetrahedron.
        volumes6 = np.linalg.det(vertices[:,1:]-vertices[:,:1])  # Evaluate the oriented volume determinant.
        if not np.all(volumes6 > 0):  # Reject a collapsed or reversed element.
            raise ValueError("The reference mesh contains an inverted or collapsed tetrahedron.")  # Stop before exporting an invalid mesh.
    used, local = np.unique(all_tets.ravel(),return_inverse=True)  # Discard construction points outside the connectivity.
    return dict(X=coordinates[used-1].T/1000, T=(local.reshape(-1,4)+1).astype(np.uint32).T,  # Return SI coordinates and one-based connectivity.
                elemCell=np.concatenate(cells), roles=np.asarray(list(roles.values())))  # Keep cell labels for the MATLAB assembly.
if __name__ == "__main__":  # Run only when invoked by the MATLAB wrapper.
    gmsh.initialize()  # Start a local meshing session.
    try:  # Guarantee cleanup after CAD or mesh errors.
        with open(sys.argv[1]) as file:  # Read the device exported by MATLAB.
            savemat(sys.argv[2], build(json.load(file)))  # Write the verified mesh.
    finally:  # Release the mesher even if export fails.
        gmsh.finalize()  # Close the native Gmsh session.
