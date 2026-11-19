def input_data(problem=2, max_it=None):
    """
    Defines input parameters for the 2D Stokes Flow Boundary Element Method (BEM)
    simulation of deformable drops in planar and constricted microchannels using
    continuous quadratic boundary elements.

    Available problems:
        1 : Straight planar microchannel (L = 10H, H = 1) with 1 drop (dad_1.geo / dad_1.msh)
        2 : Constricted microcapillary channel (L = 3H, H/h = 2) with 1 drop (dad_2.geo / dad_2.msh)
        3 : Constricted microcapillary channel (L = 6H, H/h = 4) with 1 drop (dad_3.geo / dad_3.msh)

    To change the number of time steps, edit `max_it` in the corresponding problem
    block below (or pass `max_it=<int>` when calling `input_data`).

    Boundary condition types in `bound_cond`:
        * type_x, type_y = 0 : Prescribed velocity in global x or y direction (u_x, u_y)
        * type_x, type_y = 1 : Prescribed traction in global x or y direction (t_x, t_y)
        * type_x = 2         : Prescribed parabolic Poiseuille velocity profile normal to a
                               straight inlet segment (value_x = mean velocity V_m);
                               decomposed into global (u_x, u_y) components in preproc.py
        * type_x, type_y = 3 : Fluid-fluid drop interface (two-phase capillary jump condition)

    Element boundary flags (`elem_type` stored in `conn(:, 5)` for Fortran `Processor`):
        *  0 : Solid channel wall
        * -1 : Channel inlet section (used to integrate Q_in and P_in)
        * -2 : Channel outlet section (used to integrate Q_out and P_out)
        *  1 : Deformable drop interface
    """
    name_map = {
        'dad_1': 1, 'dad1': 1, '1': 1,
        'dad_2': 2, 'dad2': 2, '2': 2,
        'dad_3': 3, 'dad3': 3, '3': 3,
    }
    if isinstance(problem, str):
        key = problem.strip().lower()
        if key in name_map:
            selected_problem = name_map[key]
        else:
            raise ValueError(f"Unknown problem selection: '{problem}'. Choose 1 ('dad_1'), 2 ('dad_2'), or 3 ('dad_3').")
    else:
        selected_problem = int(problem)

    # Common physical boundary condition mapping for channel walls, inlet, outlet, and drop
    bound_cond = {
        'wall':   {'type_x': 0, 'value_x': 0.0, 'type_y': 0, 'value_y': 0.0, 'elem_type': 0},
        'outlet': {'type_x': 1, 'value_x': 0.0, 'type_y': 1, 'value_y': 0.0, 'elem_type': -2},
        'inlet':  {'type_x': 2, 'value_x': 1.0, 'type_y': 0, 'value_y': 0.0, 'elem_type': -1},
        'drop_1': {'type_x': 3, 'value_x': 0.0, 'type_y': 3, 'value_y': 0.0, 'elem_type': 1},
    }

    n_steps_override = max_it

    if selected_problem == 1:
        # Problem 1 (dad_1): Straight plane microchannel (L = 10H, H = 1) with 1 drop (a/H = 0.50)
        file_name = 'dad_1'
        title = 'STRAIGHT MICROCHANNEL WITH DEFORMABLE DROP (dad_1)'
        dt = 0.001               # Passo de tempo (Delta t)
        max_pos = 9.5            # Posição x máxima de parada da gota
        max_it = 50              # Número de passos de tempo (use 5000 para travessia completa)
        lamb = 1.0               # Razão de viscosidade (lambda = mu_gota / mu_fluido)
        Ca = 0.1                 # Número capilar (Ca = mu * V_m / sigma)
        aH = 0.50                # Razão de confinamento (a/H = 2*R / H)

    elif selected_problem == 2:
        # Problem 2 (dad_2): Constricted microcapillary channel (L = 3H, H/h = 2) with 1 drop (a/H = 0.60)
        file_name = 'dad_2'
        title = 'CONSTRICTED MICROCHANNEL H/h=2 WITH DEFORMABLE DROP (dad_2)'
        dt = 0.001               # Passo de tempo (Delta t)
        max_pos = 3.8            # Posição x máxima de parada da gota
        max_it = 30              # Número de passos de tempo (use 8000 para travessia completa)
        lamb = 10.0              # Razão de viscosidade (lambda = mu_gota / mu_fluido)
        Ca = 1.0                 # Número capilar (Ca = mu * V_m / sigma)
        aH = 0.60                # Razão de confinamento (a/H = 2*R / H)

    elif selected_problem == 3:
        # Problem 3 (dad_3): Constricted microcapillary channel (L = 6H, H/h = 4) with 1 drop (a/H = 0.30)
        file_name = 'dad_3'
        title = 'CONSTRICTED MICROCHANNEL H/h=4 WITH DEFORMABLE DROP (dad_3)'
        dt = 0.0025              # Passo de tempo (Delta t)
        max_pos = 5.8            # Posição x máxima de parada da gota
        max_it = 50              # Número de passos de tempo (use 5000 para travessia completa)
        lamb = 10.0              # Razão de viscosidade (lambda = mu_gota / mu_fluido)
        Ca = 1.0                 # Número capilar (Ca = mu * V_m / sigma)
        aH = 0.30                # Razão de confinamento (a/H = 2*R / H)

    else:
        raise ValueError(f"Unknown problem selection: {selected_problem}. Choose 1, 2, or 3.")

    if n_steps_override is not None:
        max_it = int(n_steps_override)

    return {
        'problem': selected_problem,
        'title': title,
        'bound_cond': bound_cond,
        'file_name': file_name,
        'dt': dt,
        'max_pos': max_pos,
        'max_it': max_it,
        'lambda': lamb,
        'Ca': Ca,
        'aH': aH,
    }
