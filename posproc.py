import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

from preproc import input_data


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

    return {
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


def read_solver_outputs(results_dir='Processor'):
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

    drop_pattern = os.path.join(results_dir, 'drop_t=*.dat')
    drop_paths = sorted(
        glob.glob(drop_pattern),
        key=lambda p: int(os.path.splitext(os.path.basename(p))[0].split('=')[1]),
    )
    drop_snapshots = []
    for dp in drop_paths:
        t_ms = int(os.path.splitext(os.path.basename(dp))[0].split('=')[1])
        dxy = np.atleast_2d(np.loadtxt(dp))
        drop_snapshots.append((t_ms / 1000.0, t_ms, dxy))

    return {
        'pressure': pressure,
        'flow_error': flow_error,
        'area_error': area_error,
        'wall_contour': wall_contour,
        'drop_snapshots': drop_snapshots,
        'results_dir': results_dir,
    }


def print_validation_summary(inp, out, label='2D STOKES DROP BEM'):
    """
    Prints a formatted summary of the simulation metrics and compares conservation
    errors against the reference thresholds from the graduation thesis (UnB, 2016)
    and International Journal of Multiphase Flow (2018).
    """
    press = out['pressure']
    flow = out['flow_error']
    area = out['area_error']
    snaps = out['drop_snapshots']

    max_area_err = float(np.max(np.abs(area[:, 1])))
    mean_flow_err = float(np.mean(flow[:, 1]))

    print("=" * 78)
    print(f" {label} — NUMERICAL VALIDATION SUMMARY")
    print("=" * 78)
    print(f"Mesh Size                        : {inp['nele']} quadratic elements ({inp['nnode']} nodes)")
    print(f"Dimensionless Parameters         : lambda = {inp['lambda']:.2f}, Ca = {inp['Ca']:.2f}, a/H = {inp['aH']:.2f}")
    print(f"Simulated Time Steps             : {len(press)} steps (dt = {inp['dt']}, t_end = {press[-1, 0]:.4f})")
    print(f"Driving Pressure Delta P (start) : {press[0, 1]:.4f}")
    print(f"Driving Pressure Delta P (end)   : {press[-1, 1]:.4f}")
    print(f"Mean Flow Rate Error e_Q         : {mean_flow_err:.4f}% (reference threshold: |e_Q| <= 1.50%)")
    print(f"Max Absolute Drop Area Error |e_A|: {max_area_err:.6f}% (reference threshold: |e_A| <  0.50%)")
    if snaps:
        d0 = snaps[0][2]
        d1 = snaps[-1][2]
        print(
            f"Drop x-span evolution            : [{d0[:, 0].min():.4f}, {d0[:, 0].max():.4f}] (t={snaps[0][0]:.4f}) -> "
            f"[{d1[:, 0].min():.4f}, {d1[:, 0].max():.4f}] (t={snaps[-1][0]:.4f})"
        )
    print("=" * 78)


def plot_drop_evolution(inp, out, save_path=os.path.join('examples', 'fig1_evolucao_gota_canal.png')):
    """
    Plots the constricted microchannel contour and the Lagrangian evolution of the
    deformable drop interface, along with a zoomed-in view of the drop profiles.
    """
    n_wall = int(inp['drop_nodes'][0, 1])
    wall_xy = out['wall_contour'][:n_wall]
    snaps = out['drop_snapshots']
    cmap = plt.cm.plasma(np.linspace(0.1, 0.88, max(1, len(snaps))))

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.5, 4.4), dpi=150, gridspec_kw={'width_ratios': [1.85, 1.0]}
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
        ax2.plot(dxy[:, 0], dxy[:, 1], '-o', color=c, lw=1.6, ms=2.8, label=lbl)

    ax1.set_title(
        f'Microcapilar Constrito ($a/H={inp["aH"]:.2f}$, $\\lambda={inp["lambda"]:.1f}$, $\\mathrm{{Ca}}={inp["Ca"]:.2f}$) e Posição da Gota',
        fontsize=10.5,
        fontweight='bold',
    )
    ax1.set_xlabel('$x / H$')
    ax1.set_ylabel('$y / H$')
    ax1.set_aspect('equal', adjustable='box')
    ax1.grid(True, linestyle=':', alpha=0.45)
    ax1.legend(loc='upper right', fontsize=8, ncol=2)

    if snaps:
        t_first, t_last = snaps[0][0], snaps[-1][0]
        ax2.set_title(f'Zoom na Interface da Gota ($t = {t_first:.3f} \\to {t_last:.3f}$)', fontsize=10, fontweight='bold')
    else:
        ax2.set_title('Zoom na Interface da Gota', fontsize=10, fontweight='bold')
    ax2.set_xlabel('$x / H$')
    ax2.set_ylabel('$y / H$')
    ax2.set_aspect('equal', adjustable='box')
    ax2.grid(True, linestyle=':', alpha=0.45)
    ax2.legend(loc='upper left', fontsize=8)

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Drop evolution plot saved to '{save_path}'.")


def plot_driving_pressure(out, compare_out=None, save_path=os.path.join('examples', 'fig2_pressao_motriz.png')):
    """
    Plots the temporal evolution of the dimensionless driving pressure difference
    Delta P(t) = P_in - P_out.
    """
    press = out['pressure']

    if compare_out is not None:
        press_cmp = compare_out['pressure']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.0), dpi=150)
        ax1.plot(press[:, 0], press[:, 1], '-o', color='#1f78b4', lw=2.0, ms=3.5)
        ax1.set_title('Diferença de Pressão Motriz $\\Delta P(t)$ — Caso dad_3 ($H/h=4$)', fontsize=10.5, fontweight='bold')
        ax1.set_xlabel('Tempo adimensional $t$')
        ax1.set_ylabel('$\\Delta P = P_{\\mathrm{in}} - P_{\\mathrm{out}}$')
        ax1.grid(True, linestyle=':', alpha=0.5)

        ax2.plot(press_cmp[:, 0], press_cmp[:, 1], '-s', color='#e31a1c', lw=2.0, ms=3.5)
        ax2.set_title('Diferença de Pressão Motriz $\\Delta P(t)$ — Caso dad_2 ($H/h=2$)', fontsize=10.5, fontweight='bold')
        ax2.set_xlabel('Tempo adimensional $t$')
        ax2.set_ylabel('$\\Delta P = P_{\\mathrm{in}} - P_{\\mathrm{out}}$')
        ax2.grid(True, linestyle=':', alpha=0.5)
    else:
        fig, ax1 = plt.subplots(figsize=(7.5, 4.2), dpi=150)
        ax1.plot(press[:, 0], press[:, 1], '-o', color='#1f78b4', lw=2.0, ms=3.5)
        ax1.set_title('Diferença de Pressão Motriz $\\Delta P(t) = P_{\\mathrm{in}} - P_{\\mathrm{out}}$', fontsize=11, fontweight='bold')
        ax1.set_xlabel('Tempo adimensional $t$')
        ax1.set_ylabel('$\\Delta P = P_{\\mathrm{in}} - P_{\\mathrm{out}}$')
        ax1.grid(True, linestyle=':', alpha=0.5)

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Driving pressure plot saved to '{save_path}'.")


def plot_conservation_errors(out, compare_out=None, save_path=os.path.join('examples', 'fig3_erros_vazao_area.png')):
    """
    Plots the relative mass/flow rate conservation error e_Q(t) and the relative
    drop area conservation error e_A(t) across the simulation.
    """
    flow = out['flow_error']
    area = out['area_error']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.0), dpi=150)
    ax1.plot(flow[:, 0], flow[:, 1], '-o', color='#33a02c', lw=1.8, ms=3.2, label='Caso dad_3 (1640 nós)')
    ax1.axhline(-1.5, color='k', linestyle='--', lw=1.2, label='Limiar de convergência do TCC ($-1.5\\%$)')
    ax1.set_title('Erro Relativo de Conservação de Vazão $(Q_{\\mathrm{out}}-Q_{\\mathrm{in}})/Q_{\\mathrm{in}}$', fontsize=10.5, fontweight='bold')
    ax1.set_xlabel('Tempo adimensional $t$')
    ax1.set_ylabel('Erro de vazão (%)')
    ax1.grid(True, linestyle=':', alpha=0.5)
    ax1.legend(fontsize=8.5)

    ax2.plot(area[:, 0], area[:, 1], '-^', color='#6a3d9a', lw=1.8, ms=3.5, label='Caso dad_3 ($\\Delta t = 0.0025$)')
    if compare_out is not None:
        area_cmp = compare_out['area_error']
        ax2.plot(area_cmp[:, 0], area_cmp[:, 1], '-s', color='#ff7f00', lw=1.8, ms=3.2, label='Caso dad_2 ($\\Delta t = 0.0010$)')
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
    compare_inp=None,
    compare_out=None,
    save_path=os.path.join('examples', 'fig4_dashboard_resumo.png'),
):
    """
    Generates a 2x2 integrated summary dashboard illustrating drop trajectory/deformation,
    driving pressure evolution, and mass/area conservation errors.
    """
    n_wall = int(inp['drop_nodes'][0, 1])
    wall_xy = out['wall_contour'][:n_wall]
    snaps = out['drop_snapshots']
    press = out['pressure']
    flow = out['flow_error']
    area = out['area_error']
    cmap = plt.cm.plasma(np.linspace(0.1, 0.88, max(1, len(snaps))))

    fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.2), dpi=150)

    # Panel (a): Primary channel drop evolution
    axs[0, 0].plot(wall_xy[:, 0], wall_xy[:, 1], 'k-', lw=1.8)
    for c, (_, _, dxy) in zip(cmap, snaps):
        axs[0, 0].plot(dxy[:, 0], dxy[:, 1], '-', color=c, lw=1.6)
    axs[0, 0].set_title('(a) Evolução da Gota no Microcapilar Constrito (dad_3)', fontsize=10.5, fontweight='bold')
    axs[0, 0].set_xlabel('$x / H$')
    axs[0, 0].set_ylabel('$y / H$')
    axs[0, 0].set_aspect('equal', adjustable='box')
    axs[0, 0].grid(True, linestyle=':', alpha=0.45)

    # Panel (b): Secondary comparison case (dad_2) or zoomed drop view
    if compare_inp is not None and compare_out is not None:
        n_wall_cmp = int(compare_inp['drop_nodes'][0, 1])
        wall_xy_cmp = compare_out['wall_contour'][:n_wall_cmp]
        snaps_cmp = compare_out['drop_snapshots']
        cmap2 = plt.cm.viridis(np.linspace(0.15, 0.88, max(1, len(snaps_cmp))))
        axs[0, 1].plot(wall_xy_cmp[:, 0], wall_xy_cmp[:, 1], 'k-', lw=1.8)
        for c, (t_val, _, dxy) in zip(cmap2, snaps_cmp):
            axs[0, 1].plot(dxy[:, 0], dxy[:, 1], '-', color=c, lw=1.7, label=f'$t={t_val:.3f}$')
        axs[0, 1].set_title(f'(b) Evolução Inicial da Gota $a/H={compare_inp["aH"]:.1f}$ (dad_2 — {len(compare_out["pressure"])} passos)', fontsize=10.5, fontweight='bold')
        axs[0, 1].set_xlabel('$x / H$')
        axs[0, 1].set_ylabel('$y / H$')
        axs[0, 1].set_aspect('equal', adjustable='box')
        axs[0, 1].grid(True, linestyle=':', alpha=0.45)
        axs[0, 1].legend(loc='upper right', fontsize=8)
    else:
        for c, (t_val, _, dxy) in zip(cmap, snaps):
            axs[0, 1].plot(dxy[:, 0], dxy[:, 1], '-o', color=c, lw=1.6, ms=2.5, label=f'$t={t_val:.3f}$')
        axs[0, 1].set_title('(b) Zoom na Deformação da Interface da Gota', fontsize=10.5, fontweight='bold')
        axs[0, 1].set_xlabel('$x / H$')
        axs[0, 1].set_ylabel('$y / H$')
        axs[0, 1].set_aspect('equal', adjustable='box')
        axs[0, 1].grid(True, linestyle=':', alpha=0.45)
        axs[0, 1].legend(loc='upper right', fontsize=8)

    # Panel (c): Driving pressure evolution
    axs[1, 0].plot(press[:, 0], press[:, 1], '-o', color='#1f78b4', lw=1.8, ms=3.2)
    axs[1, 0].set_title('(c) Evolução da Queda de Pressão $\\Delta P(t)$ (dad_3)', fontsize=10.5, fontweight='bold')
    axs[1, 0].set_xlabel('Tempo adimensional $t$')
    axs[1, 0].set_ylabel('$\\Delta P$')
    axs[1, 0].grid(True, linestyle=':', alpha=0.45)

    # Panel (d): Dual-axis area and flow conservation errors
    ax_r = axs[1, 1].twinx()
    l1 = axs[1, 1].plot(
        area[:, 0], area[:, 1] * 1e3, '-^', color='#6a3d9a', lw=1.8, ms=3.2,
        label='Erro de Área ($\\times 10^{-3}\\%$)',
    )
    l2 = ax_r.plot(
        flow[:, 0], flow[:, 1], '-s', color='#33a02c', lw=1.6, ms=3.0,
        label='Erro de Vazão (%)',
    )
    axs[1, 1].set_title('(d) Conservação de Área da Gota e Vazão (dad_3)', fontsize=10.5, fontweight='bold')
    axs[1, 1].set_xlabel('Tempo adimensional $t$')
    axs[1, 1].set_ylabel('Erro de Área ($\\times 10^{-3}\\%$)', color='#6a3d9a')
    ax_r.set_ylabel('Erro de Vazão (%)', color='#33a02c')
    axs[1, 1].grid(True, linestyle=':', alpha=0.45)
    lines_all = l1 + l2
    axs[1, 1].legend(lines_all, [l.get_label() for l in lines_all], loc='center right', fontsize=8.5)

    fig.suptitle(
        'Painel Resumo Integrado 2x2 — Escoamento BEM 2D de Gotas em Microcapilar Constrito',
        fontsize=12,
        fontweight='bold',
    )
    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)
    print(f"Summary dashboard plot saved to '{save_path}'.")


def run_postprocessor(
    results_dir='Processor',
    compare_dir=os.path.join('results', 'dad_2_demo'),
    examples_dir='examples',
):
    """
    Executes the full post-processing pipeline:
      1. Reads solver input (`output.txt`) and numerical output `.dat` files;
      2. Prints validation summary tables;
      3. Generates all 4 Matplotlib figures in `examples_dir`.
    """
    inp_main = read_solver_input(os.path.join(results_dir, 'output.txt'))
    out_main = read_solver_outputs(results_dir)
    print_validation_summary(inp_main, out_main, label=f"PRIMARY CASE ({results_dir})")

    inp_cmp = None
    out_cmp = None
    if compare_dir and os.path.exists(os.path.join(compare_dir, 'pressure.dat')):
        inp_cmp = read_solver_input(os.path.join(compare_dir, 'output.txt'))
        out_cmp = read_solver_outputs(compare_dir)
        print_validation_summary(inp_cmp, out_cmp, label=f"COMPARISON CASE ({compare_dir})")

    fig1_path = os.path.join(examples_dir, 'fig1_evolucao_gota_canal.png')
    fig2_path = os.path.join(examples_dir, 'fig2_pressao_motriz.png')
    fig3_path = os.path.join(examples_dir, 'fig3_erros_vazao_area.png')
    fig4_path = os.path.join(examples_dir, 'fig4_dashboard_resumo.png')

    plot_drop_evolution(inp_main, out_main, save_path=fig1_path)
    plot_driving_pressure(out_main, compare_out=out_cmp, save_path=fig2_path)
    plot_conservation_errors(out_main, compare_out=out_cmp, save_path=fig3_path)
    plot_summary_dashboard(inp_main, out_main, compare_inp=inp_cmp, compare_out=out_cmp, save_path=fig4_path)

    return {
        'inp_main': inp_main,
        'out_main': out_main,
        'inp_cmp': inp_cmp,
        'out_cmp': out_cmp,
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
        '-c', '--compare-dir',
        default=os.path.join('results', 'dad_2_demo'),
        help='Optional secondary results directory for comparison plots (default: results/dad_2_demo).',
    )
    parser.add_argument(
        '-e', '--examples-dir',
        default='examples',
        help='Output directory for generated PNG figures (default: examples).',
    )
    args = parser.parse_args()
    run_postprocessor(
        results_dir=args.results_dir,
        compare_dir=args.compare_dir,
        examples_dir=args.examples_dir,
    )
