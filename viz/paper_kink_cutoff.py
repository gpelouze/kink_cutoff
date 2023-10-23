#!/usr/bin/env python3

import warnings

import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import papy.plot
import pyPLUTOplus as ppp

import analysis_tools as at
import cutoff_frequency
from amplitude_altitude_cut_3D_multi import AmplAltDataCollection
from phase_altitude_3D_multi import PhaseAmplitudeDataCollection


if __name__ == '__main__':

    plot_dir = 'data/plots_paper_kink_cutoff'
    period_colors = ['#525174', '#348AA7', '#5DD39E', '#BCE784']
    model_colors = ['#004488', '#bb5566', '#ddaa33', '#56b4e9']
    relax_colors = ['#000000', '#ea3d36', '#6295cb']

    fp = mpl.font_manager.FontProperties(family=['FreeSerif', 'serif'])
    plt.rcParams['figure.constrained_layout.use'] = False
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'FreeSerif'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'serif'
    plt.rcParams['mathtext.it'] = 'serif:italic'
    plt.rcParams['mathtext.bf'] = 'serif:bold'
    plt.rcParams['legend.fancybox'] = False
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['axes.formatter.limits'] = (-5, 5)

    # Data ====================================================================

    ds = ppp.PlutoDataset(
        'data/paper_kink_cutoff/104_kco_vrw_0.9995/',
        data_dir='./output',
        ns_values=[-1],
        )
    ds.flip_coordinate('x2', trans=ds.x2r.ptp())
    ds_init_2D = ppp.PlutoDataset(
        ds.ini_dir,
        data_dir=ds.data_dir,
        ns_values=[0],
        )
    ds_init_2D.flip_coordinate('x2', trans=ds_init_2D.x2r.ptp())

    wq = cutoff_frequency.WaveQuantities(ds, geometry='vertical')
    models = [
        cutoff_frequency.CutoffSpr(wq),
        cutoff_frequency.CutoffBS(wq),
        cutoff_frequency.CutoffLN(wq),
        # cutoff_frequency.CutoffPetr(wq),
        ]
    models_LN_offset = [
        cutoff_frequency.CutoffLN(wq, z0=z0)
        for z0 in u.Quantity(np.linspace(0, 2000, 4), 'km')
        ]

    wq_init_2D = cutoff_frequency.WaveQuantities(ds_init_2D, geometry='vertical')
    models_init_2D = [
        cutoff_frequency.CutoffSpr(wq_init_2D),
        cutoff_frequency.CutoffBS(wq_init_2D),
        cutoff_frequency.CutoffLN(wq_init_2D),
        # cutoff_frequency.CutoffPetr(wq_init_2D),
        ]
    models_LN_offset_init_2D = [
        cutoff_frequency.CutoffLN(wq_init_2D, z0=z0)
        for z0 in u.Quantity(np.linspace(0, 2000, 4), 'km')
        ]

    val_wq = cutoff_frequency.VALWaveQuantities()
    val_models = [
        cutoff_frequency.CutoffSpr(val_wq),
        cutoff_frequency.CutoffBS(val_wq),
        cutoff_frequency.CutoffLN(val_wq),
        # cutoff_frequency.CutoffPetr(val_wq),
        ]
    val_models_LN_offset = [
        cutoff_frequency.CutoffLN(val_wq, z0=z0)
        for z0 in u.Quantity(np.linspace(0, 2000, 4), 'km')
        ]

    class TimeDistData:
        ini_dirs = [
            'data/paper_kink_cutoff/53_kco_vrw__init_7_104__P0=200s_rerun_44',
            'data/paper_kink_cutoff/54_kco_vrw__init_7_104__P0=335s_rerun_41',
            'data/paper_kink_cutoff/51_kco_vrw__init_7_104__P0=700s_continue_42',
            'data/paper_kink_cutoff/56_kco_vrw__init_7_104__P0=2000s_rerun_43',
            ]
        dcs = []
        vx1s = []
        for ini_dir in ini_dirs:
            dc = at.CutZDataset(ini_dir, 'cut_center_pp', data_dir='./output')
            dc.flip_coordinate(trans=dc.grid.x3.xr.max())
            dc.remove_spurious_points()
            dc.list['t'] -= dc.list['t'].min()
            dc.driver_p = int(float(dc.ini['Parameters'].get('driver_p')))
            dcs.append(dc)
            vx1 = dc['vx1'] * u.Q(dc.units.velocity, 'cm s-1').to('km s-1')
            vx1s.append(vx1)

    class AmplAltData:
        filenames = [
            'data/paper_kink_cutoff/53_kco_vrw__init_7_104__P0=200s_rerun_44/plot/ampl_alt_data.npz',
            'data/paper_kink_cutoff/54_kco_vrw__init_7_104__P0=335s_rerun_41/plot/ampl_alt_data.npz',
            'data/paper_kink_cutoff/51_kco_vrw__init_7_104__P0=700s_continue_42/plot/ampl_alt_data.npz',
            'data/paper_kink_cutoff/56_kco_vrw__init_7_104__P0=2000s_rerun_43/plot/ampl_alt_data.npz',
            ]
        data = AmplAltDataCollection(filenames)

    class PhaseAltData:
        filenames = [
            'data/paper_kink_cutoff/53_kco_vrw__init_7_104__P0=200s_rerun_44/plot/phase_amplitude_data.npz',
            'data/paper_kink_cutoff/54_kco_vrw__init_7_104__P0=335s_rerun_41/plot/phase_amplitude_data.npz',
            'data/paper_kink_cutoff/51_kco_vrw__init_7_104__P0=700s_continue_42/plot/phase_amplitude_data.npz',
            'data/paper_kink_cutoff/56_kco_vrw__init_7_104__P0=2000s_rerun_43/plot/phase_amplitude_data.npz',
            ]
        data = PhaseAmplitudeDataCollection(filenames)

    # Plots ===================================================================

    class Relax2DPlot:
        u_ne = u.Quantity(ds.units.density, 'g cm-3').to('kg m-3').value
        u_ne_cross = u.Quantity(ds.units.density, 'g cm-3').to('1e-12 kg m-3').value
        u_T = u.Quantity(ds.units.temperature, 'K').to('MK').value
        u_B = u.Quantity(ds.units.magnetic_field, 'G').value

        def set_ax_color(ax, color, side):
            if side == 'right':
                ax.spines['right'].set_color(color)
                ax.spines['left'].set_visible(False)
            elif side == 'left':
                ax.spines['left'].set_color(color)
                ax.spines['right'].set_visible(False)
            else:
                raise ValueError(f'unknown side: {side}')
            ax.yaxis.label.set_color(color)
            ax.tick_params(axis='y', which='both', colors=color)

        sl_along_in = (0, np.argmin(np.abs(ds.x1 - 0)), slice(None))
        sl_along_out = (0, np.argmin(np.abs(ds.x1 - 8)), slice(None))
        sl_across = (0, slice(None), np.argmin(np.abs(ds.x2 - 30)))
        x_in = ds.x1[sl_along_in[1]]
        x_out = ds.x1[sl_along_out[1]]
        z_across = ds.x2[sl_across[2]]

        names = ['relax_strat', 'relax_init_strat']

        for ds_, name in zip([ds, ds_init_2D], names):
            plt.figure(1, clear=True)
            # setup axes
            ax_T = plt.gca()
            ax_rho = plt.twinx(ax_T)
            ax_B = plt.twinx(ax_T)
            set_ax_color(ax_T, relax_colors[0], 'left')
            set_ax_color(ax_rho, relax_colors[1], 'right')
            set_ax_color(ax_B, relax_colors[2], 'right')
            # plot
            ax_rho.plot(ds_.x2, ds_['rho'][sl_along_in] * u_ne, '-', color=relax_colors[1])
            ax_rho.plot(ds_.x2, ds_['rho'][sl_along_out] * u_ne, '--', color=relax_colors[1])
            ax_T.plot(ds_.x2, ds_['T'][sl_along_in] * u_T, '-', color=relax_colors[0])
            ax_T.plot(ds_.x2, ds_['T'][sl_along_out] * u_T, '--', color=relax_colors[0])
            ax_B.plot(ds_.x2, np.sqrt(ds_['Bx1']**2 + ds_['Bx2']**2)[sl_along_in] * u_B, '-', color=relax_colors[2])
            ax_B.plot(ds_.x2, np.sqrt(ds_['Bx1']**2 + ds_['Bx2']**2)[sl_along_out] * u_B, '--', color=relax_colors[2])
            # title
            if ds_.t[0] == 0:
                plt.title(
                    '(a) Field-aligned hydrostatic equilibrium',
                    loc='left',
                    fontsize=16.5,
                    )
            else:
                plt.title(
                    '(b) 2D magnetohydrodynamic relaxation',
                    loc='left',
                    fontsize=16.175,
                    )
            # velocity rewrite layer mark
            if ds_.t[0] > 0:
                plt.axvline(50, ls=':', lw=1, color='k')
            # legend text
            if ds_.t[0] == 0:
                ax_T.text(
                    ax_T.lines[0].get_xdata()[-1] - 5,
                    ax_T.lines[0].get_ydata()[-1],
                    '$T_\\mathrm{int}$',
                    va='bottom', ha='right',
                    color=relax_colors[0],
                    )
                ax_T.text(
                    ax_T.lines[1].get_xdata()[-1] - 5,
                    ax_T.lines[1].get_ydata()[-1],
                    '$T_\\mathrm{ext}$',
                    va='bottom', ha='right',
                    color=relax_colors[0],
                    )
                ax_rho.text(
                    ax_rho.lines[0].get_xdata()[-1] - 5,
                    ax_rho.lines[0].get_ydata()[-1] * 1.2,
                    '$\\rho_\\mathrm{int}$',
                    va='top', ha='right',
                    color=relax_colors[1],
                    )
                ax_rho.text(
                    ax_rho.lines[1].get_xdata()[-1] - 5,
                    ax_rho.lines[1].get_ydata()[-1],
                    '$\\rho_\\mathrm{ext}$',
                    va='top', ha='right',
                    color=relax_colors[1],
                    )
            else:
                ax_T.text(
                    ax_T.lines[0].get_xdata()[-1] - 5,
                    ax_T.lines[0].get_ydata()[-1],
                    '$T_\\mathrm{int}$',
                    va='top', ha='right',
                    color=relax_colors[0],
                    )
                ax_T.text(
                    ax_T.lines[1].get_xdata()[-1] - 5,
                    ax_T.lines[1].get_ydata()[-1] * 1.02,
                    '$T_\\mathrm{ext}$',
                    va='bottom', ha='right',
                    color=relax_colors[0],
                    )
                ax_rho.text(
                    ax_rho.lines[0].get_xdata()[-1] - 5,
                    ax_rho.lines[0].get_ydata()[-1] * 3.6,
                    '$\\rho_\\mathrm{int}$',
                    va='bottom', ha='right',
                    color=relax_colors[1],
                    )
                ax_rho.text(
                    ax_rho.lines[1].get_xdata()[-1] - 5,
                    ax_rho.lines[1].get_ydata()[-1] * 1.2,
                    '$\\rho_\\mathrm{ext}$',
                    va='top', ha='right',
                    color=relax_colors[1],
                    )
            ax_B.text(
                ax_B.lines[0].get_xdata()[-1] - 5,
                ax_B.lines[0].get_ydata()[-1] - .05,
                '$B_\\mathrm{int}$',
                va='top', ha='right',
                color=relax_colors[2],
                )
            ax_B.text(
                ax_B.lines[1].get_xdata()[-1] - 5,
                ax_B.lines[1].get_ydata()[-1] + .05,
                '$B_\\mathrm{ext}$',
                va='bottom', ha='right',
                color=relax_colors[2],
                )
            # scales
            plt.xscale('log')
            ax_rho.set_yscale('log')
            ax_T.set_yscale('log')
            # limits
            plt.xlim(1e-1, 1e2)
            ax_rho.set_ylim(1e-13, 5e-8)
            ax_T.set_ylim(5e-3, 5)
            # ax_B.set_ylim(10, 45)
            if ds_.t[0] == 0:
                ax_B.set_ylim(39, 44)
            else:
                ax_B.set_ylim(9, 14)
            # labels
            plt.xticks(
                ticks=[0.1, 1, 10, 100],
                labels=['0.1', '1', '10', '100'],
                )
            ax_T.set_xlabel('Altitude [Mm]')
            ax_rho.set_ylabel('Density [kg m⁻³]')
            ax_T.set_ylabel('Temperature [MK]')
            ax_B.set_ylabel('Magnetic field [G]')
            ax_B.spines['right'].set_position(('outward', 60))
            plt.tight_layout(pad=0.3)
            plt.savefig(f'{plot_dir}/{name}_vertical.pdf')

            plt.figure(1, clear=True)
            ax_T = plt.gca()
            ax_rho = plt.twinx(ax_T)
            set_ax_color(ax_rho, relax_colors[1], 'right')
            set_ax_color(ax_T, relax_colors[0], 'left')
            ax_rho.plot(ds_.x1, ds_['rho'][sl_across] * u_ne_cross, '-', color=relax_colors[1])
            ax_T.plot(ds_.x1, ds_['T'][sl_across] * u_T, '-', color=relax_colors[0], label=f'$z$ = {z_across:.0f} Mm')
            plt.xlim(0, 10)
            ax_T.set_xlabel('$x$ [Mm]')
            ax_rho.set_ylabel('Density [10⁻¹² kg m⁻³]')
            ax_T.set_ylabel('Temperature [MK]')
            ax_T.legend(loc='center right')
            plt.tight_layout(pad=0.3)
            plt.savefig(f'{plot_dir}/{name}_horizontal.pdf')

    class TimeDistPlot:
        # Prepare figure
        fig_width = 12.94 # ApJ text width
        fig_height = fig_width * .3
        figsize = (fig_width, fig_height)
        fig = plt.figure(8, clear=True, figsize=figsize)
        gs = mpl.gridspec.GridSpec(
            1, 4,
            wspace=0.1,
            left=.04, right=.99,
            top=.92, bottom=.12,
            )
        wq.dz = wq.z[1] - wq.z[0]
        wq.t_ck = np.cumsum(wq.dz / wq.ck)
        for i, (dc, vx1) in enumerate(zip(TimeDistData.dcs, TimeDistData.vx1s)):
            # Prepare data
            ix3max = dc.x3.size // 2
            vx1 = vx1[:, :ix3max]
            z = dc.x3[:ix3max]
            t = dc.t_cgs
            vx1 = vx1.to('km s-1').value
            t = t.to('s').value
            z = z.to('Mm').value
            # subplot
            ax = plt.subplot(gs[:, i])
            # Plot image
            m = papy.plot.plot_map(
                ax,
                vx1.T,
                coordinates=[t, z],
                # aspect=(t.ptp() / z.ptp()),
                aspect='auto',
                cmap='PRGn',
                norm=papy.plot.SymmetricNormalize(),
                )
            plt.plot(
                wq.t_ck.to('s').value,
                wq.z.to('Mm'),
                color='k',
                linestyle='--',
                )
            plt.ylim(0, 50)
            # plot colormap
            cb = papy.plot.colorbar(
                ax, m,
                'bottom', size='10%', pad=0.5,
                # extend=get_extend(value, **kwargs),
                extendrect=True, extendfrac=.02,
                aspect=1/30,
                )
            locator = mpl.ticker.MaxNLocator(8, steps=[1, 2, 5])
            cb.ax.xaxis.set_major_locator(locator)
            # add labels
            cb.set_label('Velocity [km s⁻¹]')
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Altitude [Mm]')
            letter = 'abcd'[i]
            ax.set_title(
                f'({letter}) $P_0 = {dc.driver_p:.0f}$ s',
                loc='left',
                fontsize=16.17,
                )
            ax.tick_params(axis='x')
            ax.tick_params(axis='y')
            cb.ax.tick_params(axis='x')
            if i != 0:
                ax.set_ylabel('')
                ax.set_yticklabels([])
        plt.savefig(f'{plot_dir}/time_dist_v.pdf')

    class AmplAltPlot:
        plt.figure(1, clear=True)
        for d, c in zip(AmplAltData.data, period_colors):
            label = f'$P_0 = {d.driver_p:.0f}$ s'
            plt.plot(
                d.z[:d.z.size//2],
                d.ampl_fit[::-1][:d.z.size//2],
                '-',
                color=c,
                label=label,
                )
        kw = dict(ha='center', va='bottom')
        plt.text(28, 14, '$P_0 = 200$ s', rotation=38, **kw)
        plt.text(31, 5.6, '$P_0 = 335$ s', rotation=30, **kw)
        plt.text(31, 2.8, '$P_0 = 700$ s', rotation=5, **kw)
        plt.text(30, 2.05, '$P_0 = 2000$ s', rotation=0.5, **kw)
        # plt.legend()
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('Velocity amplitude [km s⁻¹]')
        plt.xscale('log')
        plt.ylim(1.5, 16)
        plt.xlim(4.5e-2, 55)
        ax = plt.gca()
        # inset
        ax_zoom = ax.inset_axes([0.1, 0.5, 0.45, 0.45])
        for d, c in zip(AmplAltData.data, period_colors):
            ax_zoom.plot(
                d.z[:d.z.size//2],
                d.ampl_fit[::-1][:d.z.size//2],
                '-',
                color=c,
                label=label,
                )
        ax_zoom.set_xscale('log')
        ax_zoom.set_ylim(2, 2.6)
        ax_zoom.set_xlim(4.5e-2, 55)
        ax_zoom.set_xticks([0.1, 1, 10])
        ax_zoom.set_xticklabels(['0.1', '1', '10'])
        ax_zoom.set_xticks(
            np.concatenate([
                np.linspace(.05, .09, 5),
                np.linspace(.1, .9, 9),
                np.linspace(1, 9, 9),
                np.linspace(10, 50, 5),
                ]),
            minor=True)
        # ax_zoom.text(50, 2.2, '$P_0 = 2000$ s', ha='right')
        # ax.indicate_inset_zoom(ax_zoom, edgecolor="black")
        # ticks
        plt.xticks(
            ticks=[0.1, 1, 10],
            labels=['0.1', '1', '10'],
            )
        ax_b = plt.twiny(ax)
        ax_b.set_xscale(ax.get_xscale())
        ax_b.set_xlim(*ax.get_xlim())
        ax_b.xaxis.tick_bottom()
        ax_b.tick_params(axis='x', color='#ff000000')
        ax_b.set_xticks([0.05, 50])
        ax_b.set_xticklabels(['0.05', '50'])
        plt.tight_layout(pad=0.4)
        plt.savefig(f'{plot_dir}/ampl_alt.pdf')

    class PhaseAltPlot:
        plt.figure(1, clear=True)
        for d, c in zip(PhaseAltData.data, period_colors):
            plt.plot(
                d.z[1:d.z.size//2], 1 / d.vp[1:d.z.size//2],
                '-',
                color=c,
                label=f'$P_0 = {d.driver_p:.0f}$ s',
                )
        plt.plot(
            wq.z.to('Mm')[:wq.z.size//2],
            1 / wq.ck.to('Mm s-1')[:wq.z.size//2],
            color='C1',
            label='$1 / c_k$',
            )
        kw = dict(ha='left')
        plt.text(0.16, 58, '$1/c_k$', **kw)
        plt.text(0.15, 15, '$1/v_p$ ($P_0 = 200$ s)', rotation=-20, **kw)
        plt.text(0.15, 7, '$1/v_p$ ($P_0 = 335$ s)', rotation=-10, **kw)
        plt.text(2, 1.23, '$1/v_p$ ($P_0 = 700$ s)', rotation=-4, **kw)
        plt.text(1.8, 0.29, '$1/v_p$ ($P_0 = 2000$ s)', rotation=-14, **kw)
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$1 / v$ [s Mm⁻¹]')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1e-1, 55)
        plt.ylim(.1, 100)
        ax = plt.gca()
        # ticks
        plt.xticks(
            ticks=[0.1, 1, 10],
            labels=['0.1', '1', '10'],
            )
        ax_b = plt.twiny(ax)
        ax_b.set_xscale(ax.get_xscale())
        ax_b.set_xlim(*ax.get_xlim())
        ax_b.xaxis.tick_bottom()
        ax_b.tick_params(axis='x', color='#ff000000')
        ax_b.set_xticks([50])
        ax_b.set_xticklabels(['50'])
        plt.tight_layout(pad=0.4)
        plt.savefig(f'{plot_dir}/phase_alt.pdf')

    class CutoffFreqPlot:
        plt.figure(1, clear=True)
        ax = plt.gca()
        for m, c in zip(models, model_colors):
            if m.shortname == 'LopinNagorny2017':
                for j, m in enumerate(models_LN_offset):
                    z0_alt = m.z0
                    label = m.abbr + f' ($z_0 =$ {z0_alt.to("km"):.0f})'
                    ls = ['-', '--', '-.', ':'][j]
                    plt.plot(
                        m.z.to('Mm'),
                        m.ωc.to('s-1'),
                        color=c,
                        ls=ls,
                        label=label,
                        )
            else:
                plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=c, label=m.abbr)
        rv_thresholds = [0.2, 0.3, 0.4, 0.5]
        for i, rv_threshold in enumerate(rv_thresholds):
            N_hres = 10000
            z = []
            ω = []
            for d, c in zip(PhaseAltData.data, period_colors):
                rv = d.ck / d.vp
                z_hres = np.logspace(
                    np.log10(d.z.min()), np.log10(d.z.max()), N_hres)
                rv_hres = np.interp(z_hres, d.z, rv)
                delta = np.abs(rv_hres - rv_threshold)
                mask_precision = (delta < delta.min() * 10)
                mask_vrw = (z_hres < PhaseAltData.data.vrw_x_max_layer)
                mask_derivate = (np.gradient(rv_hres) > 0)
                mask = mask_precision & mask_vrw & mask_derivate
                if np.any(mask):
                    i_threshold = np.argmin(delta[mask])
                    z_threshold = z_hres[mask][i_threshold]
                else:
                    i_threshold = np.argmin(delta[mask_vrw])
                    z_threshold = z_hres[i_threshold]
                    msg = (f'threshold is below ck/vp for '
                           f'P0 = {d.driver_p} s and t = {rv_threshold}')
                    warnings.warn(msg)
                ω.append(d.driver_ω)
                z.append(z_threshold)
            plt.plot(
                z,
                ω,
                color='k',
                linestyle='-',
                linewidth=0.5,
                marker=['x', '+', 'o', 'D'][i],
                fillstyle='none',
                label=f'$t_r = {rv_threshold}$',
                )
        # for d, c in zip(PhaseAltData.data, period_colors):
            # plt.axhline(
                # d.driver_ω,
                # linestyle='--',
                # color=c,
                # zorder=1,
                # linewidth=1,
                # )
            # label = f'{d.driver_p:.0f} s'
            # if d.driver_p < 201:
                # label = '$P_0 = 2\\pi/\\omega_c= $' + label
            # plt.text(
                # 53, d.driver_ω,
                # label,
                # horizontalalignment='right',
                # verticalalignment='bottom',
                # )
        # tweak legend
        plt.legend(loc='lower left', ncol=2, fontsize=10)
        h, l = ax.get_legend_handles_labels()
        n_models = len(models) + len(models_LN_offset) - 1
        h_models = h[:n_models]
        l_models = l[:n_models]
        h_sims = h[n_models:]
        l_sims = l[n_models:]
        legend = ax.legend(
            h_models,
            l_models,
            title='Models',
            loc='lower left',
            bbox_to_anchor=(-.01, -.01*4/3),
            bbox_transform=ax.transAxes,
            frameon=False,
            framealpha=1,
            fontsize=10,
            )
        ax.add_artist(legend)
        ax.legend(
            h_sims,
            l_sims,
            title='Simulations',
            loc='lower left',
            bbox_to_anchor=(0.31, -.01*4/3),
            bbox_transform=ax.transAxes,
            frameon=False,
            framealpha=1,
            fontsize=10,
            )
        # labels and co
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$\\omega_c$ [s⁻¹]')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-3, 1e-1)
        plt.xlim(1e-1, 55)
        # ticks
        plt.xticks(
            ticks=[0.1, 1, 10],
            labels=['0.1', '1', '10'],
            )
        ax_b = plt.twiny(ax)
        ax_b.set_xscale(ax.get_xscale())
        ax_b.set_xlim(*ax.get_xlim())
        ax_b.xaxis.tick_bottom()
        ax_b.tick_params(axis='x', color='#ff000000')
        ax_b.set_xticks([50])
        ax_b.set_xticklabels(['50'])
        plt.tight_layout(pad=0.4)
        plt.savefig(f'{plot_dir}/cutoff_frequency.pdf')

    class VALCutoffFreqPlot:
        plt.figure(1, clear=True)
        ax = plt.gca()
        for m, c in zip(val_models, model_colors):
            if m.shortname == 'LopinNagorny2017':
                for j, m in enumerate(val_models_LN_offset):
                    z0_alt = m.z0
                    label = m.abbr + f' ($z_0 =$ {z0_alt.to("km"):.0f})'
                    ls = ['-', '--', '-.', ':'][j]
                    plt.plot(
                        m.z.to('Mm'),
                        m.ωc.to('s-1'),
                        color=c,
                        ls=ls,
                        label=label,
                        )
            else:
                plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=c, label=m.abbr)
        rv_thresholds = [0.2, 0.3, 0.4, 0.5]
        # tweak legend
        plt.legend(loc='lower left', ncol=2, fontsize=10)
        h, l = ax.get_legend_handles_labels()
        n_models = len(val_models) + len(val_models_LN_offset) - 1
        h_models = h[:n_models]
        l_models = l[:n_models]
        h_sims = h[n_models:]
        l_sims = l[n_models:]
        legend = ax.legend(
            h_models,
            l_models,
            title='Models',
            loc='lower left',
            bbox_to_anchor=(-.01, -.01*4/3),
            bbox_transform=ax.transAxes,
            frameon=False,
            framealpha=1,
            fontsize=10,
            )
        ax.add_artist(legend)
        # labels and co
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$\\omega_c$ [s⁻¹]')
        plt.xlim(-.1, 2.2)
        plt.ylim(0, 2e-2)
        # plt.xscale('log')
        # plt.yscale('log')
        # plt.ylim(1e-3, 1e-1)
        # plt.xlim(1e-1, 3)
        # # ticks
        # plt.xticks(
            # ticks=[0.1, 1, 3],
            # labels=['0.1', '1', '3'],
            # )
        ax_b = plt.twiny(ax)
        ax_b.set_xscale(ax.get_xscale())
        ax_b.set_xlim(*ax.get_xlim())
        ax_b.xaxis.tick_bottom()
        ax_b.tick_params(axis='x', color='#ff000000')
        ax_b.set_xticks([50])
        ax_b.set_xticklabels(['50'])
        plt.tight_layout(pad=0.4)
        plt.savefig(f'{plot_dir}/cutoff_frequency_val.pdf')

    class InitCutoffFreqPlot:
        plt.figure(1, clear=True)
        ax = plt.gca()
        for m1, m2, c in zip(models_init_2D, models, model_colors):
            if m1.shortname == 'LopinNagorny2017':
                for j, (m1, m2) in enumerate(zip(models_LN_offset_init_2D, models_LN_offset)):
                    z0_alt = m1.z0
                    label = m1.abbr + f' ($z_0 =$ {z0_alt.to("km"):.0f})'
                    ls = ['-', '--', '-.', ':'][j]
                    plt.plot(m1.z.to('Mm'), m1.ωc.to('s-1'), color=c, linestyle=ls, label=label)
                    # plt.plot(m2.z.to('Mm'), m2.ωc.to('s-1'), color=c, linestyle=ls, alpha=.5)
            else:
                plt.plot(m1.z.to('Mm'), m1.ωc.to('s-1'), color=c, label=m1.abbr)
                # plt.plot(m2.z.to('Mm'), m2.ωc.to('s-1'), color=c, linewidth=.5)
        # tweak legend
        plt.legend(loc='lower left', ncol=2, fontsize=10)
        h, l = ax.get_legend_handles_labels()
        n_models = len(models) + len(models_LN_offset) - 1
        h_models = h[:n_models]
        l_models = l[:n_models]
        h_sims = h[n_models:]
        l_sims = l[n_models:]
        legend = ax.legend(
            h_models,
            l_models,
            loc='lower left',
            bbox_to_anchor=(-.01, -.01*4/3),
            bbox_transform=ax.transAxes,
            frameon=True,
            framealpha=1,
            fontsize=10,
            )
        ax.add_artist(legend)
        # labels and co
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$\\omega_c$ [s⁻¹]')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-3, 1e-1)
        plt.xlim(1e-1, 55)
        # ticks
        plt.xticks(
            ticks=[0.1, 1, 10],
            labels=['0.1', '1', '10'],
            )
        ax_b = plt.twiny(ax)
        ax_b.set_xscale(ax.get_xscale())
        ax_b.set_xlim(*ax.get_xlim())
        ax_b.xaxis.tick_bottom()
        ax_b.tick_params(axis='x', color='#ff000000')
        ax_b.set_xticks([50])
        ax_b.set_xticklabels(['50'])
        plt.tight_layout(pad=0.4)
        plt.savefig(f'{plot_dir}/cutoff_frequency_init.pdf')
