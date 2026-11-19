import argparse
import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import meshio
from input_data import input_data


def ensure_mesh_exists(file_name):
    """
    Ensures `<file_name>.msh` exists and is up-to-date with `<file_name>.geo`.
    Generates the 2D quadratic mesh (`line3` + `triangle6`) from `<file_name>.geo`
    using the Gmsh Python API or the Gmsh CLI (`gmsh -2 <file>.geo -order 2`).
    """
    msh_file = file_name + '.msh'
    geo_file = file_name + '.geo'

    if os.path.exists(msh_file) and os.path.exists(geo_file):
        if os.path.getmtime(msh_file) >= os.path.getmtime(geo_file):
            return
    elif os.path.exists(msh_file) and not os.path.exists(geo_file):
        return

    if not os.path.exists(geo_file):
        raise FileNotFoundError(f"Geometry file '{geo_file}' not found.")

    try:
        import gmsh

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(geo_file)
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.setOrder(2)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(msh_file)
        gmsh.finalize()
        print(f"Mesh file '{msh_file}' generated from '{geo_file}' via Gmsh Python API.")
        return
    except Exception:
        pass

    subprocess.run(
        ['gmsh', '-2', geo_file, '-o', msh_file, '-order', '2', '-format', 'msh22'],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    print(f"Mesh file '{msh_file}' generated from '{geo_file}' via Gmsh CLI.")


def _chain_closed_loop(line_elements, seg_tags, coordinates, start_mode='origin'):
    """
    Orders a set of 3-node quadratic `line3` elements (`[n_start, n_end, n_mid]`)
    into a contiguous directed closed loop (`elem[k].n_end == elem[k+1].n_start`).

    Parameters:
        start_mode : 'origin' starts at the node closest to (0, 0) (outer channel contour);
                     'drop_front' starts at the node with maximum x-coordinate (theta = 0).
    """
    n_el = len(line_elements)
    if n_el == 0:
        return np.empty((0, 3), dtype=np.int64), np.empty((0,), dtype=np.int64)

    start_nodes = line_elements[:, 0]
    start_coords = coordinates[start_nodes, :2]

    if start_mode == 'origin':
        dists = np.hypot(start_coords[:, 0], start_coords[:, 1])
        first_idx = int(np.argmin(dists))
    elif start_mode == 'drop_front':
        # Rightmost node (theta = 0 on the drop contour)
        first_idx = int(np.argmax(start_coords[:, 0]))
    else:
        first_idx = 0

    # Map from start_node_id -> element index
    start_map = {int(line_elements[i, 0]): i for i in range(n_el)}

    ordered_indices = [first_idx]
    used = {first_idx}
    curr_end = int(line_elements[first_idx, 1])

    while len(ordered_indices) < n_el:
        next_idx = start_map.get(curr_end, None)
        if next_idx is not None and next_idx not in used:
            ordered_indices.append(next_idx)
            used.add(next_idx)
            curr_end = int(line_elements[next_idx, 1])
        else:
            # Fallback: nearest unused start node in case of floating-point duplicate vertices
            curr_pt = coordinates[curr_end, :2]
            best_j = None
            best_d = float('inf')
            for j in range(n_el):
                if j not in used:
                    d = float(np.linalg.norm(coordinates[int(line_elements[j, 0]), :2] - curr_pt))
                    if d < best_d:
                        best_d = d
                        best_j = j
            ordered_indices.append(best_j)
            used.add(best_j)
            curr_end = int(line_elements[best_j, 1])

    ordered_Idx = np.array(ordered_indices, dtype=np.int64)
    return line_elements[ordered_Idx], seg_tags[ordered_Idx]


def compute_inodes(file_name, bound_cond):
    """
    Reads a Gmsh `.msh` file and extracts nodes, quadratic boundary elements,
    1-based BEM boundary connectivity (`conec`), domain triangles, internal nodes,
    drop node ranges, and boundary condition mappings.
    """
    ensure_mesh_exists(file_name)
    mesh = meshio.read(file_name + '.msh')
    coordinates = mesh.points
    raw_line_elements = mesh.cells_dict['line3']
    triangle_elements = mesh.cells_dict['triangle6']
    raw_segments = mesh.cell_data_dict['gmsh:physical']['line3']
    edges = mesh.field_data

    bc_info = {}
    drop_keys = []
    outer_tags = set()

    for key, value in bound_cond.items():
        if key not in mesh.field_data:
            continue
        segment = int(mesh.field_data[key][0])
        elem_type = int(value.get('elem_type', 0))
        bc_info[key] = {
            'type_x': value['type_x'],
            'value_x': value['value_x'],
            'segment': segment,
            'type_y': value['type_y'],
            'value_y': value['value_y'],
            'elem_type': elem_type,
            'edge': edges[key],
        }
        if elem_type == 1 or value['type_x'] == 3:
            drop_keys.append(key)
        else:
            outer_tags.add(segment)

    # 1. Extract and chain outer channel contour elements in CCW order starting at (0, 0)
    outer_mask = np.isin(raw_segments, list(outer_tags))
    outer_elems, outer_segs = _chain_closed_loop(
        raw_line_elements[outer_mask],
        raw_segments[outer_mask],
        coordinates,
        start_mode='origin',
    )

    # 2. Extract and chain each drop contour in CW order starting at theta = 0 (front tip)
    drop_elems_list = []
    drop_segs_list = []
    for dkey in sorted(drop_keys):
        dtag = bc_info[dkey]['segment']
        dmask = (raw_segments == dtag)
        if np.any(dmask):
            d_elems, d_segs = _chain_closed_loop(
                raw_line_elements[dmask],
                raw_segments[dmask],
                coordinates,
                start_mode='drop_front',
            )
            drop_elems_list.append(d_elems)
            drop_segs_list.append(d_segs)

    if drop_elems_list:
        line_elements = np.vstack([outer_elems] + drop_elems_list)
        segments = np.concatenate([outer_segs] + drop_segs_list)
    else:
        line_elements = outer_elems
        segments = outer_segs

    # 3. Build ordered list of boundary nodes and 1-based BEM connectivity matrix
    # Note: Gmsh line3 format is [node_start, node_end, node_mid].
    # BEM local node numbering is [node_start (xi=-1), node_mid (xi=0), node_end (xi=+1)].
    ordered_boundary_nodes = []
    gmsh_to_bem = {}
    for elem in line_elements:
        for gmsh_node in (int(elem[0]), int(elem[2]), int(elem[1])):
            if gmsh_node not in gmsh_to_bem:
                ordered_boundary_nodes.append(gmsh_node)
                gmsh_to_bem[gmsh_node] = len(ordered_boundary_nodes)

    tag_to_elem_type = {info['segment']: info['elem_type'] for info in bc_info.values()}

    conec = np.zeros((line_elements.shape[0], 3), dtype=np.int64)
    elem_types = np.zeros((line_elements.shape[0],), dtype=np.int64)
    for i, elem in enumerate(line_elements):
        conec[i, 0] = gmsh_to_bem[int(elem[0])]  # Node 1 (xi = -1)
        conec[i, 1] = gmsh_to_bem[int(elem[2])]  # Node 2 (xi =  0, mid-node)
        conec[i, 2] = gmsh_to_bem[int(elem[1])]  # Node 3 (xi = +1)
        elem_types[i] = tag_to_elem_type.get(int(segments[i]), 0)

    # 4. Identify first and last BEM node of each drop contour
    drop_first_last_node = []
    elem_offset = len(outer_elems)
    for idx_drop, d_elems in enumerate(drop_elems_list, start=1):
        n_de = len(d_elems)
        first_node = int(conec[elem_offset, 0])
        last_node = int(conec[elem_offset + n_de - 1, 1])  # Midpoint of last element (highest node ID)
        drop_first_last_node.append([idx_drop, first_node, last_node])
        elem_offset += n_de

    boundary_nodes = np.unique(line_elements)
    all_nodes = np.unique(triangle_elements)
    internal_nodes = [node for node in all_nodes if node not in boundary_nodes]

    return {
        'nodes': boundary_nodes,
        'ordered_boundary_nodes': ordered_boundary_nodes,
        'conec': conec,
        'elem_types': elem_types,
        'drop_first_last_node': np.array(drop_first_last_node, dtype=np.int64),
        'internal_nodes': internal_nodes,
        'all_nodes': all_nodes,
        'coordinates': coordinates,
        'elements': line_elements,
        'triangles': triangle_elements,
        'segments': segments,
        'bc_info': bc_info,
    }


def mount_bcs(segments, bc_info, elements=None, coordinates=None):
    """
    Constructs boundary condition type (`type_bcs`) and value (`val_bcs`) arrays
    for each 3-node quadratic boundary element (6 DOFs per element).

    Supports:
      * type_x, type_y in (0, 1): Prescribed velocity (0) or traction (1)
        in global (x, y) coordinates.
      * type_x == 2: Prescribed parabolic Poiseuille velocity profile (value_x = V_m)
        normal to a straight inlet segment:
            u_n(d) = (3/2) * V_m * (1 - (2*d / H)^2)
        decomposed into global (u_x, u_y) velocity components (type = 0 in x and y).
      * type_x == 3: Deformable drop interface (assigned type = 1, value = 0.0).
    """
    num_elements = len(segments)
    type_bcs = np.zeros((num_elements, 6), dtype=np.int64)
    val_bcs = np.zeros((num_elements, 6), dtype=np.float64)

    # Precompute straight-segment geometry for any parabolic Poiseuille segment (type_x == 2)
    parab_geom = {}
    for _, bc_data in bc_info.items():
        if bc_data['type_x'] == 2:
            if elements is None or coordinates is None:
                raise ValueError("elements and coordinates are required in mount_bcs when type_x == 2.")
            seg_tag = bc_data['segment']
            seg_elems = elements[segments == seg_tag]
            if len(seg_elems) > 0:
                p_start = coordinates[int(seg_elems[0, 0]), :2]
                p_end = coordinates[int(seg_elems[-1, 1]), :2]
                delta = p_end - p_start
                length = float(np.hypot(delta[0], delta[1]))
                tg = delta / length
                # Inward normal pointing into the channel domain (CCW contour)
                norm_in = np.array([-tg[1], tg[0]], dtype=np.float64)
                center = 0.5 * (p_start + p_end)
                parab_geom[seg_tag] = (center, length, norm_in)

    for i, segment_index in enumerate(segments):
        for _, bc_data in bc_info.items():
            if bc_data['segment'] == segment_index:
                tp_x = bc_data['type_x']
                vl_x = bc_data['value_x']
                tp_y = bc_data['type_y']
                vl_y = bc_data['value_y']

                if tp_x == 2:
                    center, length, norm_in = parab_geom[int(segment_index)]
                    # Local BEM node ordering: node1 (xi=-1), node2 (xi=0), node3 (xi=+1)
                    n1, n3, n2 = int(elements[i, 0]), int(elements[i, 1]), int(elements[i, 2])
                    local_nodes = (n1, n2, n3)

                    for j, gnode in enumerate(local_nodes):
                        pt = coordinates[gnode, :2]
                        d = float(np.hypot(pt[0] - center[0], pt[1] - center[1]))
                        vel = 1.5 * vl_x * (1.0 - (2.0 * d / length) ** 2)
                        type_bcs[i, 2 * j] = 0
                        type_bcs[i, 2 * j + 1] = 0
                        val_bcs[i, 2 * j] = vel * norm_in[0]
                        val_bcs[i, 2 * j + 1] = vel * norm_in[1]

                elif tp_x == 3:
                    # Drop interface elements
                    type_bcs[i, :] = np.array([1, 1, 1, 1, 1, 1], dtype=np.int64)
                    val_bcs[i, :] = np.zeros(6, dtype=np.float64)

                else:
                    type_bcs[i, :] = np.array([tp_x, tp_y, tp_x, tp_y, tp_x, tp_y], dtype=np.int64)
                    val_bcs[i, :] = np.array([vl_x, vl_y, vl_x, vl_y, vl_x, vl_y], dtype=np.float64)
                break

    return type_bcs, val_bcs


def comp_xyplot(elem, nodes):
    """
    Interpolates coordinates along each quadratic boundary element for plotting.
    """
    npp = 8
    eta = np.linspace(-1.0, 1.0, npp)
    num_elements = elem.shape[0]

    xplot = np.zeros((num_elements, npp))
    yplot = np.zeros((num_elements, npp))

    for i in range(num_elements):
        node1, node3, node2 = elem[i]
        x1, y1 = nodes[node1][0:2]
        x2, y2 = nodes[node2][0:2]
        x3, y3 = nodes[node3][0:2]

        for k in range(npp):
            N1 = 0.5 * eta[k] * (eta[k] - 1.0)
            N2 = 1.0 - eta[k] ** 2
            N3 = 0.5 * eta[k] * (eta[k] + 1.0)
            xplot[i, k] = N1 * x1 + N2 * x2 + N3 * x3
            yplot[i, k] = N1 * y1 + N2 * y2 + N3 * y3

    return xplot, yplot


def show_geometry(coord, xplot, yplot, elem, tri6, type_bcs, val_bcs, elem_types=None, save_path='geometry.png'):
    """
    Visualizes the quadratic Gmsh mesh and global (x, y) boundary conditions of the
    2D Stokes drop problem, saving the figure to `save_path`.
    """
    num_elements = elem.shape[0]
    b_nodes = np.zeros((3 * num_elements, 2), dtype=np.float64)
    bcs = np.zeros((3 * num_elements, 4), dtype=np.float64)
    node_etype = np.zeros((3 * num_elements,), dtype=np.int64)

    for i in range(num_elements):
        local_gmsh_nodes = [elem[i, 0], elem[i, 2], elem[i, 1]]
        et = int(elem_types[i]) if elem_types is not None else 0
        for j, gnode in enumerate(local_gmsh_nodes):
            idx = 3 * i + j
            b_nodes[idx, 0] = coord[gnode, 0]
            b_nodes[idx, 1] = coord[gnode, 1]
            bcs[idx, 0] = type_bcs[i, 2 * j]
            bcs[idx, 1] = val_bcs[i, 2 * j]
            bcs[idx, 2] = type_bcs[i, 2 * j + 1]
            bcs[idx, 3] = val_bcs[i, 2 * j + 1]
            node_etype[idx] = et

    non_drop = (node_etype != 1)
    indknown_ux = np.where((bcs[:, 0] == 0) & non_drop)[0]
    indknown_tx = np.where((bcs[:, 0] == 1) & non_drop)[0]
    indknown_uy = np.where((bcs[:, 2] == 0) & non_drop)[0]
    indknown_ty = np.where((bcs[:, 2] == 1) & non_drop)[0]
    ind_drop = np.where(node_etype == 1)[0]

    x_min, x_max = np.min(coord[:, 0]), np.max(coord[:, 0])
    y_min, y_max = np.min(coord[:, 1]), np.max(coord[:, 1])
    delta_x, delta_y = x_max - x_min, y_max - y_min
    max_delta = max(delta_x, delta_y) if max(delta_x, delta_y) > 0 else 1.0
    max_bc_val = np.max(np.abs(bcs[:, [1, 3]]))
    arrow_scale = (8.0 * max_bc_val / max_delta) if max_bc_val > 0 else 10.0

    aspect = max(1.6, min(3.2, delta_x / max(delta_y, 1e-6)))
    fig, ax = plt.subplots(figsize=(3.8 * aspect, 4.2))
    ax.grid(True, linestyle='--', alpha=0.4)

    # Plot the underlying triangular Gmsh mesh
    if tri6 is not None and len(tri6) > 0:
        tri3 = tri6[:, :3]
        ax.triplot(coord[:, 0], coord[:, 1], tri3, color="lightgray", linewidth=0.4)

    # Plot quadratic boundary element curves
    for i in range(num_elements):
        et = int(elem_types[i]) if elem_types is not None else 0
        color = '#d95f02' if et == 1 else 'black'
        lw = 1.6 if et == 1 else 1.0
        line_segments = LineCollection(
            [list(zip(xplot[i, :], yplot[i, :]))],
            linestyles='solid',
            colors=color,
            linewidths=lw,
        )
        ax.add_collection(line_segments)

    # Plot global x-direction velocity vectors (red)
    if len(indknown_ux) > 0:
        nonzero_ux = indknown_ux[np.abs(bcs[indknown_ux, 1]) > 1e-12]
        if len(nonzero_ux) > 0:
            ax.quiver(
                b_nodes[nonzero_ux, 0], b_nodes[nonzero_ux, 1],
                bcs[nonzero_ux, 1], np.zeros_like(bcs[nonzero_ux, 1]),
                color="red", angles='xy', scale_units='xy', scale=arrow_scale,
                width=0.003,
            )

    # Plot markers for BCs and drop nodes
    if len(indknown_tx) > 0:
        ax.plot(b_nodes[indknown_tx, 0], b_nodes[indknown_tx, 1], 'bx', markersize=5, label='Known $t_x, t_y$ (Outlet)')
    if len(indknown_ux) > 0:
        ax.plot(b_nodes[indknown_ux, 0], b_nodes[indknown_ux, 1], 'r.', markersize=4, label='Known $u_x, u_y$ (Wall / Inlet)')
    if len(ind_drop) > 0:
        ax.plot(b_nodes[ind_drop, 0], b_nodes[ind_drop, 1], 'o', color='#d95f02', markersize=3.5, label='Drop Interface ($\\Gamma_d$)')

    ax.set_aspect("equal")
    ax.set_xlim(x_min - 0.12 * max_delta, x_max + 0.08 * max_delta)
    ax.set_ylim(y_min - 0.12 * delta_y, y_max + 0.18 * delta_y)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Gmsh Quadratic Mesh and Boundary Conditions (2D Stokes BEM)')
    ax.legend(loc='upper right', fontsize=8.5, framealpha=0.9)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Geometry plot saved to '{save_path}'.")
    plt.close(fig)


def create_input_file(filename, computed_data, type_bcs, val_bcs, input_params):
    """
    Writes the formatted input file (`Processor/output.txt`) for the Fortran
    2D Stokes BEM solver (`Gota.exe`), matching the exact format expected by
    `Processor/read_input.f90`.
    """
    nodes = computed_data['coordinates']
    ordered_boundary_nodes = computed_data['ordered_boundary_nodes']
    conec = computed_data['conec']
    elem_types = computed_data['elem_types']
    drop_first_last_node = computed_data['drop_first_last_node']

    nnode = len(ordered_boundary_nodes)
    nele = conec.shape[0]
    ndrop = len(drop_first_last_node)

    dt = float(input_params['dt'])
    max_pos = float(input_params['max_pos'])
    max_it = int(input_params['max_it'])
    lamb = float(input_params['lambda'])
    Ca = float(input_params['Ca'])
    aH = float(input_params['aH'])

    os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)

    with open(filename, 'w', encoding='utf-8') as f:
        # 1. Header: nnode, nele, ndrop
        f.write(f"{nnode:10d}{nele:10d}{ndrop:10d}        #number of nodes, elements and drops\n")
        f.write(" \n")

        # 2. Time-stepping parameters: dt, max_pos, max_it
        f.write(f"{dt:18.7f}{max_pos:18.7f}{max_it:10d}   #time step, max drop position and max iteration\n")
        f.write(" \n")

        # 3. Node coordinates
        f.write(" NODES COORDINATES (i# node, x, y)\n")
        f.write(" \n")
        for bem_id, gmsh_node_id in enumerate(ordered_boundary_nodes, start=1):
            x_val = float(nodes[gmsh_node_id, 0])
            y_val = float(nodes[gmsh_node_id, 1])
            f.write(f"{float(bem_id):18.7f}{x_val:18.7f}{y_val:18.7f}\n")

        # 4. Element connectivity (1-based: elem_id, node1, node2, node3, elem_type)
        f.write(" \n")
        f.write(" CONNECTIVITY (i# elem, 1st nodes, 2nd node, 3th node, type)\n")
        f.write(" \n")
        for i in range(nele):
            f.write(
                f"{i + 1:10d}{int(conec[i, 0]):10d}{int(conec[i, 1]):10d}"
                f"{int(conec[i, 2]):10d}{int(elem_types[i]):10d}\n"
            )

        # 5. Element boundary conditions (13 columns per element)
        f.write(" \n")
        f.write(" ELEM BOUNDARY CONDITIONS (i# elem, nd type x, vlaue x, nd type y, value y)\n")
        f.write(" \n")
        for i in range(nele):
            row_vals = [
                float(i + 1),
                float(type_bcs[i, 0]), float(val_bcs[i, 0]),
                float(type_bcs[i, 1]), float(val_bcs[i, 1]),
                float(type_bcs[i, 2]), float(val_bcs[i, 2]),
                float(type_bcs[i, 3]), float(val_bcs[i, 3]),
                float(type_bcs[i, 4]), float(val_bcs[i, 4]),
                float(type_bcs[i, 5]), float(val_bcs[i, 5]),
            ]
            f.write("".join(f"{v:18.7f}" for v in row_vals) + "\n")

        # 6. Drop node ranges and physical properties
        f.write(" \n")
        f.write(" DROPS FIRST AND LAST NODES (i# drop, first node, last node)\n")
        f.write(" \n")
        for i in range(ndrop):
            f.write(
                f"{int(drop_first_last_node[i, 0]):7d}"
                f"{int(drop_first_last_node[i, 1]):7d}"
                f"{int(drop_first_last_node[i, 2]):7d}\n"
            )

        f.write(" \n")
        f.write(" DROPS PROPERTIES (Lambda, Ca, aH)\n")
        f.write(" \n")
        f.write(f"{lamb:18.7f}{Ca:18.7f}{aH:18.7f}\n")
        f.write(" \n")

    print(f"Input file '{filename}' generated successfully ({nnode} nodes, {nele} elements, {ndrop} drop(s)).")


def run_preprocessor(problem=2, max_it=None, output_path=None, plot_path='geometry.png'):
    """
    Runs the complete Python preprocessing pipeline for the selected problem:
      1. Reads problem definition from `input_data(problem, max_it=max_it)`;
      2. Generates/loads `<file_name>.msh` from `<file_name>.geo` via Gmsh + meshio;
      3. Evaluates nodal boundary conditions (including Poiseuille parabolic inlet);
      4. Plots the Gmsh mesh and boundary conditions to `plot_path`;
      5. Writes the Fortran solver input file (`Processor/output.txt`).
    """
    inp_data = input_data(problem, max_it=max_it)
    computed_data = compute_inodes(inp_data['file_name'], inp_data['bound_cond'])
    type_bcs, val_bcs = mount_bcs(
        computed_data['segments'],
        computed_data['bc_info'],
        computed_data['elements'],
        computed_data['coordinates'],
    )
    xplot, yplot = comp_xyplot(computed_data['elements'], computed_data['coordinates'])
    show_geometry(
        computed_data['coordinates'],
        xplot,
        yplot,
        computed_data['elements'],
        computed_data['triangles'],
        type_bcs,
        val_bcs,
        elem_types=computed_data['elem_types'],
        save_path=plot_path,
    )

    if output_path is None:
        output_path = os.path.join('Processor', 'output.txt')

    create_input_file(
        output_path,
        computed_data,
        type_bcs,
        val_bcs,
        inp_data,
    )
    return inp_data, computed_data, type_bcs, val_bcs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='2D Stokes BEM Python Preprocessor (Gmsh + meshio)')
    parser.add_argument(
        '-p', '--problem',
        default=2,
        help="Problem number or name: 1 ('dad_1'), 2 ('dad_2', default), or 3 ('dad_3').",
    )
    parser.add_argument(
        '-o', '--output',
        default=os.path.join('Processor', 'output.txt'),
        help="Path to the output file read by the Fortran solver (default: Processor/output.txt).",
    )
    parser.add_argument(
        '--plot',
        default='geometry.png',
        help="Path to save the mesh and boundary conditions figure (default: geometry.png).",
    )
    parser.add_argument(
        '-n', '--max-it',
        type=int,
        default=None,
        help="Optional override for number of time steps (default: uses max_it from input_data.py).",
    )
    args = parser.parse_args()
    run_preprocessor(problem=args.problem, max_it=args.max_it, output_path=args.output, plot_path=args.plot)
