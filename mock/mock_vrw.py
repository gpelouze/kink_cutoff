#!/usr/bin/env python3

import argparse

import etframes
import matplotlib.pyplot as plt
import numpy as np
import papy


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    # pluto parameters
    p.add_argument('tstop', type=float)
    p.add_argument('--dt1', type=float, default=0)
    p.add_argument('--dt2', type=float, default=0)
    p.add_argument('--dt3', type=float, default=0)
    p.add_argument('--dt4', type=float, default=0)
    p.add_argument('--dt5', type=float, default=0)
    p.add_argument('--dt6', type=float, default=0)
    p.add_argument('--avfmin', type=float, default=0.9)
    p.add_argument('--avlmin', type=float, default=0.8)
    p.add_argument('--xmaxl', type=float, default=50.)
    p.add_argument('--L', type=float, default=100)
    # mock parameters
    p.add_argument('--Nt', type=int, default=500)
    p.add_argument('--Nx', type=int, default=200)
    p.add_argument('--plot-symbolic', action='store_true')
    p.add_argument('--real-scale', action='store_true')
    args = p.parse_args()

    all_time = np.linspace(0, args.tstop, args.Nt)
    x = np.linspace(0, args.L, args.Nx)

    dt1 = args.dt1
    dt2 = args.dt2
    dt3 = args.dt3
    dt4 = args.dt4
    dt5 = args.dt5
    dt6 = args.dt6
    av_fulldom_min = args.avfmin
    av_layer_min = args.avlmin
    x_max_layer = args.xmaxl

    av_fulldom_min_print = av_fulldom_min
    av_layer_min_print = av_layer_min
    if not args.real_scale:
        av_fulldom_min = np.clip(av_fulldom_min, None, 0.9)
        av_layer_min = np.clip(av_layer_min, None, 0.8)
    if av_fulldom_min_print < av_layer_min_print:
        av_fulldom_min, av_layer_min = av_layer_min, av_fulldom_min

    # mock

    av = np.full((args.Nt, args.Nx), np.nan)
    av_fulldom = av.copy()
    av_layer = av.copy()
    bv_layer = av.copy()

    t1 = dt1
    t2 = t1 + dt2
    t3 = t2 + dt3
    t4 = t3 + dt4
    t5 = t4 + dt5
    t6 = t5 + dt6

    for it, g_time in enumerate(all_time):

        for ix, _ in enumerate(x):

            # alpha_fulldom
            if g_time < t1:
                av_fulldom[it, ix] = 1.
            elif g_time < t2:
                av_fulldom[it, ix] = 1. - (g_time - t1) * (1. - av_fulldom_min) / dt2
            elif g_time < t3:
                av_fulldom[it, ix] = av_fulldom_min
            elif g_time < t4:
                av_fulldom[it, ix] = av_fulldom_min + (g_time - t3) * (1. - av_fulldom_min) / dt4
            elif g_time < t5:
                av_fulldom[it, ix] = 1.
            elif g_time < t6:
                if (dt4 == 0) and (dt5 == 0):
                    av_fulldom[it, ix] = av_fulldom_min - (g_time - t5) * av_fulldom_min / dt6
                else:
                    av_fulldom[it, ix] = 1. - (g_time - t5) / dt6
            else:
                av_fulldom[it, ix] = 0.

            # alpha_layer
            if g_time < t5:
                av_layer[it, ix] = 0.
            elif g_time < t6:
                av_layer[it, ix] = (g_time - t5) / dt6
            else:
                av_layer[it, ix] = 1.

            # beta_layer
            if x[ix] < x_max_layer:
                bv_layer[it, ix] = (
                    av_layer_min +
                    (1 - av_layer_min) / x_max_layer * x[ix]
                    )
            else:
                bv_layer[it, ix] = 1.

    av = av_fulldom + av_layer*bv_layer
    assert np.all(av <= 1.)

    # plot result

    dt = {i+1: v for i, v in enumerate([dt1, dt2, dt3, dt4, dt5, dt6])}
    t = {i+1: v for i, v in enumerate([0, t1, t2, t3, t4, t5, t6])}
    dt_nonzero = {i: dt[i] for i in dt.keys() if dt[i] > 0}
    t_nonzero = {i: t[i+1] for i in dt.keys() if dt[i] > 0}

    t_ticks = {}
    if all_time.min() not in t_nonzero.values():
        t_ticks[0] = all_time.min()
    t_ticks.update(t_nonzero)
    if all_time.max() not in t_nonzero.values():
        t_ticks['end'] = all_time.max()
    if args.plot_symbolic:
        t_ticklabels = {i: f'$t_\\mathrm{{{i}}}$' for i, t in t_ticks.items()}
        if 0 in t_ticklabels:
            t_ticklabels[0] = 0
    else:
        t_ticklabels = {i: f'{t:g}' for i, t in t_ticks.items()}
    t_ticks = list(t_ticks.values())
    t_ticklabels = list(t_ticklabels.values())

    x_ticks = [0, x_max_layer, args.L]
    if args.plot_symbolic:
        x_ticklabels = ['0 (apex)', '$x_{\\mathrm{max},l}$', 'L']
    else:
        x_ticklabels = [0., float(x_max_layer), float(args.L)]
    x_ticks_all = x_ticks.copy()
    x_ticklabels_all = x_ticklabels.copy()
    if t5 >= args.tstop:
        del x_ticks[1]
        del x_ticklabels[1]

    av_fulldom_plot = papy.num.almost_identical(av_fulldom, 1e-10, axis=1)
    av_layer_plot = papy.num.almost_identical(av_layer, 1e-10, axis=1)
    plt.clf()
    plt.plot(
        all_time,
        av_layer_plot,
        '--',
        color='#008B72',
        label='$\\alpha_{v,l}$',
        )
    plt.plot(
        all_time,
        av_fulldom_plot,
        '-',
        color='#ddaa33',
        label='$\\alpha_{v,f}$',
        )
    for i, t in t_nonzero.items():
        if t > 0:
            plt.plot([t, t], [0, 1], color='k', alpha=.2)
        if args.plot_symbolic:
            plt.text(t-dt[i]/2, -0.075, f'$\\mathrm{{d}}t_{i}$', ha='center', va='top', color='gray')
    plt.legend(frameon=False)
    plt.xlabel('Time')
    plt.ylabel('$\\alpha_v$')
    etframes.add_range_frame()
    plt.yticks(
        ticks=(0, av_fulldom_min, 1),
        labels=('0',
                ('$\\alpha_{v,f,\\mathrm{min}}$' if args.plot_symbolic
                 else f'{av_fulldom_min_print:g}'),
                1),
        )
    plt.xticks(ticks=t_ticks, labels=t_ticklabels)
    plt.savefig('data/mock_vrw_alpha.pdf')

    bv_layer_plot = papy.num.almost_identical(bv_layer, 1e-10, axis=0)
    plt.clf()
    plt.plot(
        x,
        bv_layer_plot,
        '--',
        color='#008B72',
        label='$\\beta_{v,l}$',
        )
    plt.plot([x_max_layer, x_max_layer], [av_layer_min, 1], color='k', alpha=.2)
    plt.xlabel('Position along the loop')
    plt.ylabel('$\\beta_v$')
    etframes.add_range_frame()
    plt.yticks(
        ticks=(av_layer_min, 1),
        labels=(
            ('$\\alpha_{v,l,\\mathrm{min}}$' if args.plot_symbolic
             else av_layer_min_print),
            1,
            ),
        )
    plt.xticks(ticks=x_ticks_all, labels=x_ticklabels_all)
    plt.savefig('data/mock_vrw_beta.pdf')

    plt.clf()
    m = papy.plot.plot_map(
        plt.gca(),
        av.T,
        coordinates=[all_time, x],
        aspect=(all_time.ptp() / x.ptp()),
        cmap='gray',
        )

    if 1 in dt_nonzero:
        plt.text(dt1/2, args.L/2,
                 '1.0',
                 ha='center', va='center')
    if 3 in dt_nonzero:
        plt.text(t3-dt3/2, args.L/2,
                 f'{av_fulldom_min_print:g}',
                 color='white',
                 ha='center', va='center')
    if 5 in dt_nonzero:
        plt.text(t5-dt5/2, args.L/2,
                 '1.0',
                 ha='center', va='center')
    if t6 < args.tstop:
        plt.text(t6+(args.tstop-t6)/2,
                 x_max_layer + (args.L-x_max_layer)/2,
                 f'1.0',
                 ha='center', va='center')
        plt.text(t6+(args.tstop-t6)/2, 0,
                 f'{av_layer_min_print:g}',
                 color='white',
                 ha='center', va='bottom')

    cb = plt.colorbar(m, label='$\\alpha_v$', pad=0)
    plt.xlabel('Time')
    plt.ylabel('Position along the loop')
    plt.xticks(ticks=t_ticks, labels=t_ticklabels)
    plt.yticks(ticks=x_ticks, labels=x_ticklabels)
    cb.set_ticks((0, av_layer_min, av_fulldom_min, 1))
    cb.set_ticklabels((
        '0',
        ('$\\alpha_{v,l,\\mathrm{min}}$' if args.plot_symbolic
         else av_layer_min_print),
        ('$\\alpha_{v,f,\\mathrm{min}}$' if args.plot_symbolic
         else av_fulldom_min_print),
        1))
    plt.savefig('data/mock_vrw_av.pdf')
