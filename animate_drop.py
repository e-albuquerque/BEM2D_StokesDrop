import argparse
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from posproc import read_solver_input, read_solver_outputs


def _configure_ffmpeg():
    """
    Configures Matplotlib to use system `ffmpeg` if available, or falls back to
    the bundled binary provided by `imageio-ffmpeg`.
    Returns True if an FFmpeg executable is ready, False otherwise.
    """
    if shutil.which('ffmpeg') is not None:
        return True
    try:
        import imageio_ffmpeg

        ffmpeg_exe = imageio_ffmpeg.get_ffmpeg_exe()
        if ffmpeg_exe and os.path.exists(ffmpeg_exe):
            plt.rcParams['animation.ffmpeg_path'] = ffmpeg_exe
            return True
    except ImportError:
        pass
    return animation.FFMpegWriter.isAvailable()


def create_drop_animation(
    results_dir='Processor',
    output_path=os.path.join('examples', 'drop_animation.mp4'),
    fps=15,
    dpi=150,
    also_gif=False,
):
    """
    Generates an animation video (MP4 and optionally GIF) showing the deformable
    drop flowing through the microcapillary channel synchronized with the driving
    pressure difference Delta P(t).
    """
    inp = read_solver_input(os.path.join(results_dir, 'output.txt'))
    out = read_solver_outputs(results_dir, inp=inp)

    snaps = out['drop_snapshots']
    if not snaps:
        raise RuntimeError(f"No 'drop_t=*.dat' snapshot files found in '{results_dir}'.")

    press = out['pressure']
    area = out['area_error']
    n_wall = int(inp['drop_nodes'][0, 1])
    wall_xy = out['wall_contour'][:n_wall]
    x_max_wall = float(np.max(wall_xy[:, 0]))

    # Compute drop centroids and synchronize pressure/area values for every snapshot
    snap_times = np.array([s[0] for s in snaps], dtype=np.float64)
    snap_press = np.interp(snap_times, press[:, 0], press[:, 1])
    snap_area_err = np.interp(snap_times, area[:, 0], area[:, 1])
    centroids = np.array([[float(np.mean(s[2][:-1, 0])), float(np.mean(s[2][:-1, 1]))] for s in snaps])

    # Set up 2-panel figure: top = channel + drop, bottom = synchronized Delta P(t)
    fig, (ax_ch, ax_dp) = plt.subplots(
        2, 1,
        figsize=(10.4, 6.2),
        dpi=dpi,
        gridspec_kw={'height_ratios': [1.55, 1.0]},
    )

    # --- Top panel: microchannel geometry ---
    ax_ch.plot(wall_xy[:, 0], wall_xy[:, 1], 'k-', lw=2.2, zorder=4)
    if x_max_wall > 5.0:
        ax_ch.fill_between([2.0, 4.0], [0.0, 0.0], [0.375, 0.375], color='#dcdcdc', zorder=1)
        ax_ch.fill_between([2.0, 4.0], [0.625, 0.625], [1.0, 1.0], color='#dcdcdc', zorder=1)
    elif x_max_wall >= 2.9:
        ax_ch.fill_between([1.0, 2.0], [0.0, 0.0], [0.25, 0.25], color='#dcdcdc', zorder=1)
        ax_ch.fill_between([1.0, 2.0], [0.75, 0.75], [1.0, 1.0], color='#dcdcdc', zorder=1)

    ax_ch.set_xlim(-0.08 * x_max_wall, 1.05 * x_max_wall)
    ax_ch.set_ylim(-0.12, 1.15)
    ax_ch.set_aspect('equal', adjustable='box')
    ax_ch.set_xlabel('$x / H$', fontsize=10)
    ax_ch.set_ylabel('$y / H$', fontsize=10)
    ax_ch.set_title(
        f'Escoamento BEM 2D de Gota Deformável ({inp["case_name"]} — '
        f'$a/H={inp["aH"]:.2f}$, $\\lambda={inp["lambda"]:.1f}$, $\\mathrm{{Ca}}={inp["Ca"]:.2f}$)',
        fontsize=11.5,
        fontweight='bold',
    )
    ax_ch.grid(True, linestyle=':', alpha=0.45)

    # Dynamic artists on top panel
    d0 = snaps[0][2]
    drop_patch, = ax_ch.fill(d0[:, 0], d0[:, 1], facecolor='#f4a582', edgecolor='none', alpha=0.75, zorder=2)
    drop_line, = ax_ch.plot(d0[:, 0], d0[:, 1], '-', color='#b2182b', lw=2.0, zorder=3)
    trail_line, = ax_ch.plot([centroids[0, 0]], [centroids[0, 1]], '--', color='#2166ac', lw=1.2, alpha=0.7, zorder=2)
    info_text = ax_ch.text(
        0.015, 0.93, '',
        transform=ax_ch.transAxes,
        fontsize=9.2,
        verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#cccccc', alpha=0.92),
        zorder=5,
    )

    # --- Bottom panel: synchronized driving pressure Delta P(t) ---
    ax_dp.plot(press[:, 0], press[:, 1], '-', color='#cccccc', lw=1.5, label='Curva completa $\\Delta P(t)$')
    dp_trace, = ax_dp.plot([], [], '-', color='#1f78b4', lw=2.2, label='Evolução atual')
    dp_dot, = ax_dp.plot([], [], 'o', color='#e31a1c', ms=6.5, zorder=5)

    dp_min, dp_max = float(np.min(press[:, 1])), float(np.max(press[:, 1]))
    dp_pad = max(0.5, 0.12 * (dp_max - dp_min))
    ax_dp.set_xlim(float(press[0, 0]), float(press[-1, 0]))
    ax_dp.set_ylim(dp_min - dp_pad, dp_max + dp_pad)
    ax_dp.set_xlabel('Tempo adimensional $t$', fontsize=10)
    ax_dp.set_ylabel('Pressão motriz $\\Delta P(t)$', fontsize=10)
    ax_dp.grid(True, linestyle=':', alpha=0.5)
    ax_dp.legend(loc='upper right', fontsize=8.5, framealpha=0.9)

    fig.tight_layout()

    def update(frame_idx):
        t_val, _, dxy = snaps[frame_idx]
        drop_patch.set_xy(dxy)
        drop_line.set_data(dxy[:, 0], dxy[:, 1])
        trail_line.set_data(centroids[:frame_idx + 1, 0], centroids[:frame_idx + 1, 1])

        mask = press[:, 0] <= t_val + 1e-9
        dp_trace.set_data(press[mask, 0], press[mask, 1])
        dp_dot.set_data([t_val], [snap_press[frame_idx]])

        info_text.set_text(
            f'$t = {t_val:.3f}$  |  '
            f'$\\Delta P = {snap_press[frame_idx]:.2f}$  |  '
            f'$x_{{\\mathrm{{ponta}}}} = {dxy[:, 0].max():.3f}$  |  '
            f'$e_A = {snap_area_err[frame_idx]:+.3f}\\%$'
        )
        return drop_patch, drop_line, trail_line, dp_trace, dp_dot, info_text

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(snaps),
        interval=1000.0 / fps,
        blit=False,
    )

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    ext = os.path.splitext(output_path)[1].lower()

    saved_files = []
    if ext == '.mp4':
        if _configure_ffmpeg():
            writer = animation.FFMpegWriter(
                fps=fps,
                codec='libx264',
                extra_args=['-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-pix_fmt', 'yuv420p'],
            )
            anim.save(output_path, writer=writer, dpi=dpi)
            print(f"MP4 animation saved to '{output_path}' ({len(snaps)} frames @ {fps} fps).")
            saved_files.append(output_path)
        else:
            fallback_gif = os.path.splitext(output_path)[0] + '.gif'
            print("[!] FFmpeg not found; falling back to GIF via PillowWriter.")
            anim.save(fallback_gif, writer=animation.PillowWriter(fps=fps), dpi=dpi)
            print(f"GIF animation saved to '{fallback_gif}' ({len(snaps)} frames @ {fps} fps).")
            saved_files.append(fallback_gif)
    elif ext == '.gif':
        anim.save(output_path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
        print(f"GIF animation saved to '{output_path}' ({len(snaps)} frames @ {fps} fps).")
        saved_files.append(output_path)
    else:
        raise ValueError(f"Unsupported output format '{ext}'. Use '.mp4' or '.gif'.")

    if also_gif and ext != '.gif':
        gif_path = os.path.splitext(output_path)[0] + '.gif'
        anim.save(gif_path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
        print(f"GIF animation saved to '{gif_path}' ({len(snaps)} frames @ {fps} fps).")
        saved_files.append(gif_path)

    plt.close(fig)
    return saved_files


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate MP4/GIF animation of the 2D BEM drop simulation')
    parser.add_argument(
        '-r', '--results-dir',
        default='Processor',
        help='Directory containing output.txt and *.dat results (default: Processor).',
    )
    parser.add_argument(
        '-o', '--output',
        default=os.path.join('examples', 'drop_animation.mp4'),
        help='Output animation path (.mp4 or .gif; default: examples/drop_animation.mp4).',
    )
    parser.add_argument(
        '--fps',
        type=int,
        default=15,
        help='Frames per second (default: 15).',
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=150,
        help='Resolution in DPI (default: 150).',
    )
    parser.add_argument(
        '--gif',
        action='store_true',
        help='Also export a GIF copy alongside the MP4 video.',
    )
    args = parser.parse_args()
    create_drop_animation(
        results_dir=args.results_dir,
        output_path=args.output,
        fps=args.fps,
        dpi=args.dpi,
        also_gif=args.gif,
    )
