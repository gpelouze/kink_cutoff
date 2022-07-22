#!/usr/bin/env python3

"""
Generate velocity rewrite layer parameter using new simplified formulation

$\alpha_v(t, z) = 1 - \beta_v(z) \gamma_v(t)$

This formulation is equivalent to the one in `mock_vrw.py` when
$\alpha_{v,f,\mathrm{min}} = 1$.
"""

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
    p.add_argument('--avmin', type=float, default=0.8)
    p.add_argument('--xmaxl', type=float, default=50.)
    p.add_argument('--L', type=float, default=100)
    # mock parameters
    p.add_argument('--Nt', type=int, default=500)
    p.add_argument('--Nx', type=int, default=200)
    p.add_argument('--split-norm', action='store_true')
    p.add_argument('--plot-symbolic', action='store_true')
    args = p.parse_args()

    all_time = np.linspace(0, args.tstop, args.Nt)
    x = np.linspace(0, args.L, args.Nx)

    dt1 = args.dt1
    dt2 = args.dt2
    av_min = args.avmin
    x_max_layer = args.xmaxl
    domain_length = args.L

    # mock

    av = np.full((args.Nt, args.Nx), np.nan)
    bv = av.copy()
    cv = av.copy()

    t1 = dt1
    t2 = t1 + dt2

    for it, g_time in enumerate(all_time):

        for ix, _ in enumerate(x):

            # beta_v,l
            if x[ix] < x_max_layer:
                bv[it, ix] = (
                    (1 - av_min) *
                    (domain_length - x_max_layer - x[ix]) /
                    (domain_length - x_max_layer)
                    )
            else:
                bv[it, ix] = 0.

            # gamma_v,l
            if g_time <= t1:
                cv[it, ix] = 0.
            elif g_time <= t2:
                cv[it, ix] = (g_time - t1) / dt2
            else:
                cv[it, ix] = 1.

    av = 1 - bv * cv
    assert np.all(av <= 1.)

    # plot result

    dt = {i+1: v for i, v in enumerate([dt1, dt2])}
    t = {i+1: v for i, v in enumerate([0, t1, t2])}
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
    if t2 >= args.tstop:
        del x_ticks[1]
        del x_ticklabels[1]

    norm_imshow = plt.matplotlib.colors.Normalize(
        vmin=av_min, vmax=1)

    cv_plot = papy.num.almost_identical(cv, 1e-10, axis=1)
    plt.figure(clear=True, constrained_layout=False)
    plt.plot(
        all_time,
        cv_plot,
        '-',
        color='#ddaa33',
        )
    for i, t in t_nonzero.items():
        if t > 0:
            plt.plot([t, t], [0, 1], color='k', alpha=.2)
        if args.plot_symbolic:
            plt.text(t-dt[i]/2, -0.075, f'$\\mathrm{{d}}t_{i}$', ha='center', va='top', color='gray')
    plt.xlabel('Time')
    plt.ylabel('$\\gamma_v$')
    etframes.add_range_frame()
    plt.yticks(
        ticks=(0, 1),
        )
    plt.xticks(ticks=t_ticks, labels=t_ticklabels)
    plt.tight_layout()
    plt.savefig('data/new_mock_vrw_gamma.pdf')

    bv_plot = papy.num.almost_identical(bv, 1e-10, axis=0)
    plt.figure(clear=True, constrained_layout=False)
    plt.plot(
        x,
        bv_plot,
        '--',
        color='#008B72',
        )
    plt.plot([x_max_layer, x_max_layer], [0, 1 - av_min], color='k', alpha=.2)
    plt.xlabel('Position along the loop')
    plt.ylabel('$\\beta_v$')
    etframes.add_range_frame()
    plt.yticks(
        ticks=(0, 1 - av_min),
        labels=(
            0,
            ('$1 - \\alpha_{v,\\mathrm{min}}$' if args.plot_symbolic
             else f'{1 - av_min:g}'),
            ),
        )
    plt.xticks(ticks=x_ticks_all, labels=x_ticklabels_all)
    plt.tight_layout()
    plt.savefig('data/new_mock_vrw_beta.pdf')

    plt.figure(clear=True, constrained_layout=False)
    m = papy.plot.plot_map(
        plt.gca(),
        av.T,
        coordinates=[all_time, x],
        aspect=(all_time.ptp() / x.ptp()),
        cmap='gray',
        norm=norm_imshow,
        )

    if 1 in dt_nonzero:
        plt.text(t1-dt1/2, args.L/2,
                 f'{1:g}',
                 bbox=dict(facecolor='white', edgecolor='black'),
                 ha='center', va='center')
    if t2 < args.tstop:
        plt.text(t2+(args.tstop-t2)/2,
                 x_max_layer + (args.L-x_max_layer)/2,
                 f'{1:g}',
                 bbox=dict(facecolor='white', edgecolor='black'),
                 ha='center', va='center')
        plt.text(t2 + (args.tstop-t2) / 2, 0,
                 f'{av_min:g}',
                 color='black',
                 bbox=dict(facecolor='white', edgecolor='black'),
                 ha='center', va='bottom')

    cb = plt.colorbar(m, label='$\\alpha_v$', pad=0.02)
    plt.xlabel('Time')
    plt.ylabel('Position along the loop')
    plt.xticks(ticks=t_ticks, labels=t_ticklabels)
    plt.yticks(ticks=x_ticks, labels=x_ticklabels)
    cb.set_ticks((av_min, 1))
    cb.set_ticklabels((
        ('$\\alpha_{v,\\mathrm{min}}$' if args.plot_symbolic
         else av_min),
        1))
    plt.tight_layout()
    plt.savefig('data/new_mock_vrw_av.pdf')
