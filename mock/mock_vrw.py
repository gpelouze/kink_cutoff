#!/usr/bin/env python3

import etframes
import matplotlib.pyplot as plt
import numpy as np
import papy


if __name__ == '__main__':

    # mock parameters
    Nt = 500
    Nx = 200

    # pluto parameters
    tstop = 1400.
    L = 100
    all_time = np.linspace(0, tstop, Nt)
    x = np.linspace(0, L, Nx)

    # pluto userdef parameters
    dt1 = 200.
    dt2 = 200.
    dt3 = 200.
    dt4 = 200.
    dt5 = 200
    dt6 = 200
    vrw_av_fulldom_min = .8
    vrw_av_layer_min = .5
    vrw_x_max_layer = 50.

    # mock
    av = np.full((Nt, Nx), np.nan)
    av_fulldom = av.copy()
    av_layer = av.copy()
    bv_layer = av.copy()

    t1 = dt1
    t2 = t1 + dt2
    t3 = t2 + dt3
    t4 = t3 + dt4
    t5 = t4 + dt5
    t6 = t5 + dt6
    av_fulldom_min = vrw_av_fulldom_min
    av_layer_min = vrw_av_layer_min
    x_max_layer = vrw_x_max_layer

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
    for i, (dt, t) in enumerate(zip([dt1, dt2, dt3, dt4, dt5, dt6], [0, t1, t2, t3, t4, t5])):
        if dt > 0:
            if t > 0:
                plt.plot([t, t], [0, 1], color='k', alpha=.2)
            plt.text(t+dt/2, -0.075, f'$\\mathrm{{d}}t_{i+1}$', ha='center', va='top', color='gray')
    plt.legend(frameon=False)
    plt.xlabel('Time')
    plt.ylabel('$\\alpha_v$')
    etframes.add_range_frame()
    plt.yticks(
        ticks=(0, av_fulldom_min, 1),
        labels=('0', '$\\alpha_{v,f,\\mathrm{min}}$', 1),
        )
    plt.xticks(
        ticks=(all_time.min(), t1, t2, t3, t4, t5, t6, all_time.max()),
        labels=('0', '$t_1$', '$t_2$', '$t_3$', '$t_4$', '$t_5$', '$t_6$', '$t_\\mathrm{stop}$'),
        )
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
    plt.plot([vrw_x_max_layer, vrw_x_max_layer], [av_layer_min, 1], color='k', alpha=.2)
    plt.xlabel('Position along the loop')
    plt.ylabel('$\\beta_v$')
    etframes.add_range_frame()
    plt.yticks(
        ticks=(av_layer_min, 1),
        labels=('$\\alpha_{v,l,\\mathrm{min}}$', 1),
        )
    plt.xticks(
        ticks=(x.min(), x_max_layer, x.max()),
        labels=('0 (apex)', '$x_{\\mathrm{max},l}$', 'L/2'),
        )
    plt.savefig('data/mock_vrw_beta.pdf')

    plt.clf()
    m = papy.plot.plot_map(
        plt.gca(),
        av.T,
        coordinates=[all_time, x],
        aspect=(all_time.ptp() / x.ptp()),
        cmap='gray',
        )
    cb = plt.colorbar(m, label='$\\alpha_v$', pad=0)
    plt.xlabel('Time')
    plt.ylabel('Position along the loop')
    plt.xticks(
        ticks=(all_time.min(), t1, t2, t3, t4, t5, t6, all_time.max()),
        labels=('0', '$t_1$', '$t_2$', '$t_3$', '$t_4$', '$t_5$', '$t_6$', '$t_\\mathrm{stop}$'),
        )
    plt.yticks(
        ticks=(x.min(), x_max_layer, x.max()),
        labels=('0', '$x_{\\mathrm{max},l}$', 'L/2'),
        )
    cb.set_ticks((0, av_layer_min, av_fulldom_min, 1))
    cb.set_ticklabels(('0', '$\\alpha_{v,l,\\mathrm{min}}$', '$\\alpha_{v,f,\\mathrm{min}}$', 1))
    plt.savefig('data/mock_vrw_av.pdf')
