#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import papy


if __name__ == '__main__':

    # mock parameters
    Nt = 500
    Nx = 200

    # pluto parameters
    tstop = 1500.
    L = 100
    all_time = np.linspace(0, tstop, Nt)
    x = np.linspace(0, L, Nx)

    # pluto userdef parameters
    dt1 = 300.
    dt2 = 100.
    dt3 = 300.
    dt4 = 100.
    dt5 = 300
    dt6 = 100
    vrw_av_fulldom_min = .9
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
        av_fulldom_plot,
        'r-',
        label='$\\alpha_\\mathrm{relax}$',
        )
    plt.plot(
        all_time,
        av_layer_plot,
        'g--',
        label='$\\alpha_\\mathrm{layer}$',
        )
    for tmark in [t1, t2, t3, t4, t5, t6]:
        plt.axvline(tmark, alpha=.2)
    plt.legend()
    plt.xlabel('Time [$t_0$]')
    plt.ylabel('$\\alpha$')
    plt.savefig('data/mock_vrw_alpha.pdf')

    bv_layer_plot = papy.num.almost_identical(bv_layer, 1e-10, axis=0)
    plt.clf()
    plt.plot(
        x,
        bv_layer_plot,
        'g--',
        label='$\\beta_\\mathrm{layer}$',
        )
    plt.axvline(vrw_x_max_layer, alpha=.2)
    plt.legend()
    plt.xlabel('x [Mm]')
    plt.ylabel('$\\beta$')
    plt.savefig('data/mock_vrw_beta.pdf')

    plt.clf()
    m = papy.plot.plot_map(
        plt.gca(),
        av.T,
        coordinates=[all_time, x],
        aspect=(all_time.ptp() / x.ptp()),
        vmin=0,
        vmax=1,
        )
    plt.colorbar(m, label='$\\alpha_v$')
    plt.xlabel('Time [$t_0$]')
    plt.ylabel('x [Mm]')
    plt.savefig('data/mock_vrw_av.pdf')
