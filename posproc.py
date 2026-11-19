import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt


def read_solver_input(filename=os.path.join('Processor', 'output.txt')):
    """
    Reads the Fortran 2D Stokes BEM solver input file (`output.txt`) and returns
    mesh dimensions, time-stepping parameters, boundary node coordinates, element
    connectivity (`conn`), boundary conditions (`bc`), drop node indices, and
    physical dimensionless parameters (`lambda`, `Ca`, `aH`).
    """
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.read().splitlines()

    nnode, nele, ndrop = map(int, lines[0].split('#')[0].split()[:3])
    s2 = lines[2].split('#')[0].split()
    dt, max_pos, max_it = float(s2[0]), float(s2[1]), int(s2[2])

    xnode = np.array(
        [[float(x) for x in lines[6 + i].split()[:3]] for i in range(nnode)],
        dtype=np.float64,
    )
    c0 = 6 + nnode + 3
    conn = np.array(
        [[int(x) for x in lines[c0 + i].split()[:5]] for i in range(nele)],
        dtype=np.int64,
    )
    b0 = c0 + nele + 3
    bc = np.array(
        [[float(x) for x in lines[b0 + i].split()[:13]] for i in range(nele)],
        dtype=np.float64,
    )
    d0 = b0 + nele + 3
    drop_nodes = np.array(
        [[int(x) for x in lines[d0 + i].split()[:3]] for i in range(ndrop)],
        dtype=np.int64,
    )
    p0 = lines[d0 + ndrop + 3].split()
    lamb, Ca, aH = float(p0[0]), float(p0[1]), float(p0[2])

    if nele == 470:
        case_name = 'dad_1'
    elif nele == 220:
        case_name = 'dad_2'
    elif nele == 820:
        case_name = 'dad_3'
    else:
        case_name = f'{nele}_elem'

    return {
        'case_name': case_name,
        'nnode': nnode,
        'nele': nele,
        'ndrop': ndrop,
        'dt': dt,
        'max_pos': max_pos,
        'max_it': max_it,
        'xnode': xnode,
        'conn': conn,
        'bc': bc,
        'drop_nodes': drop_nodes,
        'lambda': lamb,
        'Ca': Ca,
        'aH': aH,
    }


def read_solver_outputs(results_dir='Processor', inp=None):
    """
    Reads numerical results exported by the 2D Stokes BEM solver (`Gota.exe`):
      - `pressure.dat`    : Dimensionless time t and driving pressure difference Delta P(t)
      - `flow_error.dat`  : Dimensionless time t and relative outlet flow rate error e_Q(t) (%)
      - `area_error.dat`  : Dimensionless time t and relative drop area error e_A(t) (%)
      - `wall_contour.dat`: Channel wall nodal coordinates (x, y)
      - `drop_t=*.dat`    : Time snapshots of the deforming drop contour (x, y)
    """
    pressure = np.atleast_2d(np.loadtxt(os.path.join(results_dir, 'pressure.dat')))
    flow_error = np.atleast_2d(np.loadtxt(os.path.join(results_dir, 'flow_error.dat')))
    area_error = np.atleast_2d(np.loadtxt(os.path.join(results_dir, 'area_error.dat')))
    wall_contour = np.atleast_2d(np.loadtxt(os.path.join(results_dir, 'wall_contour.dat')))

    # Valid t_ms values corresponding to the active simulation in pressure.dat
    valid_t_ms = set(np.round(pressure[:, 0] * 1000.0).astype(int).tolist())
    n_expected_drop_rows = (int(inp['drop_nodes'][0, 2]) - int(inp['drop_nodes'][0, 1]) + 2) if inp is not None else None

    drop_pattern = os.path.join(results_dir, 'drop_t=*.dat')
    drop_paths = sorted(
        glob.glob(drop_pattern),
        key=lambda p: int(os.path.splitext(os.path.basename(p))[0].split('=')[1]),
    )
    drop_snapshots = []
    for dp in drop_paths:
        t_ms = int(os.path.splitext(os.path.basename(dp))[0].split('=')[1])
        if valid_t_ms and t_ms not in valid_t_ms:
            continue
        dxy = np.atleast_2d(np.loadtxt(dp))
        if n_expected_drop_rows is not None and dxy.shape[0] != n_expected_drop_rows:
            continue
        drop_snapshots.append((t_ms / 1000.0, t_ms, dxy))

    return {
        'pressure': pressure,
        'flow_error': flow_error,
        'area_error': area_error,
        'wall_contour': wall_contour,
        'drop_snapshots': drop_snapshots,
        'results_dir': results_dir,
    }


def _select_snapshots(snaps, max_count=10):
    """
    Selects up to `max_count` evenly spaced snapshots (always including the first
    and the last snapshot) so plots and legends remain clean for long runs.
    """
    n = len(snaps)
    if n <= max_count:
        return snaps
    indices = np.unique(np.round(np.linspace(0, n - 1, max_count)).astype(int))
    return [snaps[i] for i in indices]


def print_validation_summary(inp, out, label=None):
    """
    Prints a formatted summary of the simulation metrics and compares conservation
    errors against the reference thresholds from the graduation thesis (UnB, 2016)
    and International Journal of Multiphase Flow (2018).
    """
    press = out['pressure']
    flow = out['flow_error']
    area = out['area_error']
    snaps = out['drop_snapshots']

    if label is None:
        label = f"2D STOKES DROP BEM ({inp['case_name']})"

    max_area_err = float(np.max(np.abs(area[:, 1])))
    mean_flow_err = float(np.mean(flow[:, 1]))

    print("=" * 78)
    print(f" {label} — NUMERICAL VALIDATION SUMMARY")
    print("=" * 78)
    print(f"Case / Mesh Size                 : {inp['case_name']} — {inp['nele']} quadratic elements ({inp['nnode']} nodes)")
    print(f"Dimensionless Parameters         : lambda = {inp['lambda']:.2f}, Ca = {inp['Ca']:.2f}, a/H = {inp['aH']:.2f}")
    print(f"Simulated Time Steps             : {len(press)} steps (dt = {inp['dt']}, t_end = {press[-1, 0]:.4f})")
    print(f"Driving Pressure Delta P (min/max): {press[:, 1].min():.4f} / {press[:, 1].max():.4f} (start={press[0, 1]:.4f}, end={press[-1, 1]:.4f})")
    print(f"Mean Flow Rate Error e_Q         : {mean_flow_err:.4f}% (reference threshold: |e_Q| <= 1.50%)")
    print(f"Max Absolute Drop Area Error |e_A|: {max_area_err:.6f}% (reference threshold: |e_A| <  0.50%)")
    if snaps:
        d0 = snaps[0][2]
        d1 = snaps[-1][2]
        print(
            f"Drop x-span evolution            : [{d0[:, 0].min():.4f}, {d0[:, 0].max():.4f}] (t={snaps[0][0]:.3f}) -> "
            f"[{d1[:, 0].min():.4f}, {d1[:, 0].max():.4f}] (t={snaps[-1][0]:.3f})"
        )
    print("=" * 78)


def plot_drop_evolution(inp, out, save_path=os.path.join('examples', 'fig1_evolucao_gota_canal.png')):
    """
    Plots the constricted microchannel contour and the Lagrangian evolution of the
    deformable drop interface, along with a zoomed-in view of representative profiles.
    """
    n_wall = int(inp['drop_nodes'][0, 1])
    wall_xy = out['wall_contour'][:n_wall]
    snaps = _select_snapshots(out['drop_snapshots'], max_count=10)
    n_steps = len(out['pressure'])
    cmap = plt.cm.plasma(np.linspace(0.08, 0.90, max(1, len(snaps))))

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(13.0, 4.5), dpi=150, gridspec_kw={'width_ratios': [1.65, 1.15]}
    )
    ax1.plot(wall_xy[:, 0], wall_xy[:, 1], 'k-', lw=2.0, label='Parede do Microcapilar')

    x_max_wall = float(np.max(wall_xy[:, 0]))
    if x_max_wall > 5.0:
        ax1.fill_between([2.0, 4.0], [0.0, 0.0], [0.375, 0.375], color='#e0e0e0', alpha=0.7)
        ax1.fill_between([2.0, 4.0], [0.625, 0.625], [1.0, 1.0], color='#e0e0e0', alpha=0.7)
    elif x_max_wall >= 2.9:
        ax1.fill_between([1.0, 2.0], [0.0, 0.0], [0.25, 0.25], color='#e0e0e0', alpha=0.7)
        ax1.fill_between([1.0, 2.0], [0.75, 0.75], [1.0, 1.0], color='#e0e0e0', alpha=0.7)

    for c, (t_val, _, dxy) in zip(cmap, snaps):
        lbl = f'$t = {t_val:.3f}$'
        ax1.plot(dxy[:, 0], dxy[:, 1], '-', color=c, lw=1.8, label=lbl)
        ax2.plot(dxy[:, 0], dxy[:, 1], '-o', color=c, lw=1.5, ms=2.4, label=lbl)

    ax1.set_title(
        f'Evolução da Gota no Microcapilar ({inp["case_name"]} — {n_steps} passos, $a/H={inp["aH"]:.2f}$, $\\lambda={inp["lambda"]:.1f}$, $\\mathrm{{Ca}}={inp["Ca"]:.2f}$)',
        fontsize=10.2,
        fontweight='bold',
    )
    ax1.set_xlabel('$x / H$')
    ax1.set_ylabel('$y / H$')
    ax1.set_aspect('equal', adjustable='box')
    ax1.grid(True, linestyle=':', alpha=0.45)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.26), fontsize=8.0, ncol=6, framealpha=0.95)

    if snaps:
        t_first, t_last = snaps[0][0], snaps[-1][0]
        ax2.set_title(f'Perfis da Interface da Gota ($t = {t_first:.3f} \\to {t_last:.3f}$)', fontsize=10.2, fontweight='bold')
    else:
        ax2.set_title('Perfis da Interface da Gota', fontsize=10.2, fontweight='bold')
    ax2.set_xlabel('$x / H$')
    ax2.set_ylabel('$y / H$')
    ax2.set_aspect('equal', adjustable='box')
    ax2.grid(True, linestyle=':', alpha=0.45)
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.32), fontsize=8.0, ncol=5, framealpha=0.95)

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Drop evolution plot saved to '{save_path}'.")


def plot_driving_pressure(inp, out, save_path=os.path.join('examples', 'fig2_pressao_motriz.png')):
    """
    Plots the temporal evolution of the dimensionless driving pressure difference
    Delta P(t) = P_in - P_out and the extra pressure drop Delta P^+(t) = Delta P(t) - Delta P(0).
    """
    press = out['pressure']
    n_steps = len(press)
    dp_extra = press[:, 1] - press[0, 1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.2), dpi=150)
    ax1.plot(press[:, 0], press[:, 1], '-', color='#1f78b4', lw=2.0)
    ax1.set_title(
        f'Pressão Motriz $\\Delta P(t) = P_{{\\mathrm{{in}}}} - P_{{\\mathrm{{out}}}}$ ({inp["case_name"]} — {n_steps} passos)',
        fontsize=10.5,
        fontweight='bold',
    )
    ax1.set_xlabel('Tempo adimensional $t$')
    ax1.set_ylabel('$\\Delta P(t)$')
    ax1.grid(True, linestyle=':', alpha=0.5)

    ax2.plot(press[:, 0], dp_extra, '-', color='#e31a1c', lw=2.0)
    ax2.axhline(0.0, color='k', linestyle='--', lw=1.0, alpha=0.6)
    ax2.set_title(
        f'Variação de Pressão Adicional $\\Delta P^+(t) = \\Delta P(t) - \\Delta P(t_0)$ ({inp["case_name"]})',
        fontsize=10.5,
        fontweight='bold',
    )
    ax2.set_xlabel('Tempo adimensional $t$')
    ax2.set_ylabel('$\\Delta P^+(t)$')
    ax2.grid(True, linestyle=':', alpha=0.5)

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Driving pressure plot saved to '{save_path}'.")


def plot_conservation_errors(inp, out, save_path=os.path.join('examples', 'fig3_erros_vazao_area.png')):
    """
    Plots the relative mass/flow rate conservation error e_Q(t) and the relative
    drop area conservation error e_A(t) across all simulated time steps.
    """
    flow = out['flow_error']
    area = out['area_error']
    n_steps = len(flow)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.2), dpi=150)
    ax1.plot(flow[:, 0], flow[:, 1], '-', color='#33a02c', lw=1.8, label=f'Caso {inp["case_name"]} ({inp["nnode"]} nós, {n_steps} passos)')
    ax1.axhline(-1.5, color='k', linestyle='--', lw=1.2, label='Limiar de convergência do TCC ($-1.5\\%$)')
    ax1.set_title('Erro Relativo de Conservação de Vazão $(Q_{\\mathrm{out}}-Q_{\\mathrm{in}})/Q_{\\mathrm{in}}$', fontsize=10.5, fontweight='bold')
    ax1.set_xlabel('Tempo adimensional $t$')
    ax1.set_ylabel('Erro de vazão (%)')
    ax1.grid(True, linestyle=':', alpha=0.5)
    ax1.legend(fontsize=8.5)

    ax2.plot(area[:, 0], area[:, 1], '-', color='#6a3d9a', lw=1.8, label=f'Caso {inp["case_name"]} ($\\Delta t = {inp["dt"]}$, {n_steps} passos)')
    ax2.set_title('Erro Relativo de Conservação da Área da Gota $(A - A_0)/A_0$', fontsize=10.5, fontweight='bold')
    ax2.set_xlabel('Tempo adimensional $t$')
    ax2.set_ylabel('Erro de área (%)')
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.legend(fontsize=8.5)

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Conservation errors plot saved to '{save_path}'.")


def plot_summary_dashboard(
    inp,
    out,
    save_path=os.path.join('examples', 'fig4_dashboard_resumo.png'),
):
    """
    Generates a 2x2 integrated summary dashboard for the active simulation in
    `Processor/`, illustrating:
      (a) Full microchannel contour and drop trajectory across all steps;
      (b) Deformed drop interface profiles across the constriction (all steps);
      (c) Driving pressure evolution Delta P(t);
      (d) Drop area and flow rate conservation errors.
    """
    n_wall = int(inp['drop_nodes'][0, 1])
    wall_xy = out['wall_contour'][:n_wall]
    snaps_all = out['drop_snapshots']
    snaps = _select_snapshots(snaps_all, max_count=10)
    press = out['pressure']
    flow = out['flow_error']
    area = out['area_error']
    n_steps = len(press)
    cmap = plt.cm.plasma(np.linspace(0.08, 0.90, max(1, len(snaps))))

    fig, axs = plt.subplots(2, 2, figsize=(12.8, 8.2), dpi=150)

    # Panel (a): Full channel + drop evolution across all steps
    axs[0, 0].plot(wall_xy[:, 0], wall_xy[:, 1], 'k-', lw=1.8)
    x_max_wall = float(np.max(wall_xy[:, 0]))
    if x_max_wall > 5.0:
        axs[0, 0].fill_between([2.0, 4.0], [0.0, 0.0], [0.375, 0.375], color='#e0e0e0', alpha=0.7)
        axs[0, 0].fill_between([2.0, 4.0], [0.625, 0.625], [1.0, 1.0], color='#e0e0e0', alpha=0.7)
    elif x_max_wall >= 2.9:
        axs[0, 0].fill_between([1.0, 2.0], [0.0, 0.0], [0.25, 0.25], color='#e0e0e0', alpha=0.7)
        axs[0, 0].fill_between([1.0, 2.0], [0.75, 0.75], [1.0, 1.0], color='#e0e0e0', alpha=0.7)

    for c, (_, _, dxy) in zip(cmap, snaps):
        axs[0, 0].plot(dxy[:, 0], dxy[:, 1], '-', color=c, lw=1.7)
    axs[0, 0].set_title(
        f'(a) Trajetória da Gota no Microcapilar ({inp["case_name"]} — {n_steps} passos)',
        fontsize=10.5,
        fontweight='bold',
    )
    axs[0, 0].set_xlabel('$x / H$')
    axs[0, 0].set_ylabel('$y / H$')
    axs[0, 0].set_aspect('equal', adjustable='box')
    axs[0, 0].grid(True, linestyle=':', alpha=0.45)

    # Panel (b): Detailed drop deformation profiles across all simulated steps
    for c, (t_val, _, dxy) in zip(cmap, snaps):
        axs[0, 1].plot(dxy[:, 0], dxy[:, 1], '-', color=c, lw=1.8, label=f'$t={t_val:.3f}$')
    axs[0, 1].set_title(
        f'(b) Deformação da Gota $a/H={inp["aH"]:.2f}$ ({inp["case_name"]} — {n_steps} passos)',
        fontsize=10.5,
        fontweight='bold',
    )
    axs[0, 1].set_xlabel('$x / H$')
    axs[0, 1].set_ylabel('$y / H$')
    axs[0, 1].set_aspect('equal', adjustable='box')
    axs[0, 1].grid(True, linestyle=':', alpha=0.45)
    axs[0, 1].legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.34),
        fontsize=8.0,
        ncol=5,
        framealpha=0.95,
    )

    # Panel (c): Driving pressure evolution
    axs[1, 0].plot(press[:, 0], press[:, 1], '-', color='#1f78b4', lw=2.0)
    axs[1, 0].set_title(
        f'(c) Evolução da Queda de Pressão $\\Delta P(t)$ ({inp["case_name"]} — {n_steps} passos)',
        fontsize=10.5,
        fontweight='bold',
    )
    axs[1, 0].set_xlabel('Tempo adimensional $t$')
    axs[1, 0].set_ylabel('$\\Delta P(t)$')
    axs[1, 0].grid(True, linestyle=':', alpha=0.45)

    # Panel (d): Dual-axis area and flow conservation errors
    ax_r = axs[1, 1].twinx()
    l1 = axs[1, 1].plot(
        area[:, 0], area[:, 1], '-', color='#6a3d9a', lw=1.8,
        label='Erro de Área (%)',
    )
    l2 = ax_r.plot(
        flow[:, 0], flow[:, 1], '--', color='#33a02c', lw=1.6,
        label='Erro de Vazão (%)',
    )
    axs[1, 1].set_title(
        f'(d) Conservação de Área da Gota e Vazão ({inp["case_name"]} — {n_steps} passos)',
        fontsize=10.5,
        fontweight='bold',
    )
    axs[1, 1].set_xlabel('Tempo adimensional $t$')
    axs[1, 1].set_ylabel('Erro de Área (%)', color='#6a3d9a')
    ax_r.set_ylabel('Erro de Vazão (%)', color='#33a02c')
    axs[1, 1].grid(True, linestyle=':', alpha=0.45)
    lines_all = l1 + l2
    axs[1, 1].legend(lines_all, [l.get_label() for l in lines_all], loc='best', fontsize=8.5)

    fig.suptitle(
        f'Painel Resumo Integrado 2x2 — Escoamento BEM 2D de Gota ({inp["case_name"]}, {n_steps} passos, $a/H={inp["aH"]:.2f}$, $\\lambda={inp["lambda"]:.1f}$, $\\mathrm{{Ca}}={inp["Ca"]:.2f}$)',
        fontsize=11.5,
        fontweight='bold',
    )
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Summary dashboard plot saved to '{save_path}'.")


def run_postprocessor(
    results_dir='Processor',
    examples_dir='examples',
):
    """
    Executes the full post-processing pipeline for the simulation results in `results_dir`:
      1. Reads solver input (`output.txt`) and numerical output `.dat` files;
      2. Prints validation summary table;
      3. Generates all 4 Matplotlib figures in `examples_dir`.
    """
    inp_main = read_solver_input(os.path.join(results_dir, 'output.txt'))
    out_main = read_solver_outputs(results_dir, inp=inp_main)
    print_validation_summary(inp_main, out_main)

    fig1_path = os.path.join(examples_dir, 'fig1_evolucao_gota_canal.png')
    fig2_path = os.path.join(examples_dir, 'fig2_pressao_motriz.png')
    fig3_path = os.path.join(examples_dir, 'fig3_erros_vazao_area.png')
    fig4_path = os.path.join(examples_dir, 'fig4_dashboard_resumo.png')

    plot_drop_evolution(inp_main, out_main, save_path=fig1_path)
    plot_driving_pressure(inp_main, out_main, save_path=fig2_path)
    plot_conservation_errors(inp_main, out_main, save_path=fig3_path)
    plot_summary_dashboard(inp_main, out_main, save_path=fig4_path)

    return {
        'inp_main': inp_main,
        'out_main': out_main,
        'figures': [fig1_path, fig2_path, fig3_path, fig4_path],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='2D Stokes BEM Python Post-Processor (Matplotlib)')
    parser.add_argument(
        '-r', '--results-dir',
        default='Processor',
        help='Directory containing output.txt and *.dat results (default: Processor).',
    )
    parser.add_argument(
        '-e', '--examples-dir',
        default='examples',
        help='Output directory for generated PNG figures (default: examples).',
    )
    args = parser.parse_args()
    run_postprocessor(
        results_dir=args.results_dir,
        examples_dir=args.examples_dir,
    )
