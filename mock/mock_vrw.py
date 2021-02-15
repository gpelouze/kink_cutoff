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
    vrv_t_ramp = 100.
    vrv_t_prs = 300.
    vrv_t_wave = 300.
    vrv_t_restore = 300.
    vrv_av_relax = .9
    vrv_av_layer_min = .0
    vrv_x_layer_max = 50.

    # mock
    av = np.full((Nt, Nx), np.nan)
    av_relax = av.copy()
    av_layer = av.copy()
    bv_layer = av.copy()

    t0 = 0.
    t1 = t0 + vrv_t_prs
    t2 = t1 + vrv_t_ramp
    t3 = t2 + vrv_t_wave
    t4 = t3 + vrv_t_ramp
    t5 = t4 + vrv_t_restore
    t6 = t5 + vrv_t_ramp

    for it, g_time in enumerate(all_time):

        for ix, _ in enumerate(x):

            # beta_layer
            if x[ix] < vrv_x_layer_max:
                bv_layer[it, ix] = (
                    vrv_av_layer_min +
                    (1 - vrv_av_layer_min) / vrv_x_layer_max * x[ix]
                    )
            else:
                bv_layer[it, ix] = 1.

            # alpha_relax
            if g_time < t1:
                av_relax[it, ix] = 1.
            elif g_time < t2:
                av_relax[it, ix] = 1. - (g_time - t1) * (1. - vrv_av_relax) / vrv_t_ramp
            elif g_time < t3:
                av_relax[it, ix] = vrv_av_relax
            elif g_time < t4:
                av_relax[it, ix] = vrv_av_relax + (g_time - t3) * (1. - vrv_av_relax) / vrv_t_ramp
            elif g_time < t5:
                av_relax[it, ix] = 1.
            elif g_time < t6:
                av_relax[it, ix] = 1. - (g_time - t5) / vrv_t_ramp
            else:
                av_relax[it, ix] = 0.

            # alpha_layer
            if g_time < t5:
                av_layer[it, ix] = 0.
            elif g_time < t6:
                av_layer[it, ix] = (g_time - t5) / vrv_t_ramp
            else:
                av_layer[it, ix] = 1.

    av = av_relax + av_layer*bv_layer
    assert np.all(av <= 1.)

    # plot result
    av_relax_plot = papy.num.almost_identical(av_relax, 1e-10, axis=1)
    av_layer_plot = papy.num.almost_identical(av_layer, 1e-10, axis=1)
    plt.clf()
    plt.plot(
        all_time,
        av_relax_plot,
        'r-',
        label='$\\alpha_\\mathrm{relax}$',
        )
    plt.plot(
        all_time,
        av_layer_plot,
        'g--',
        label='$\\alpha_\\mathrm{layer}$',
        )
    for tmark in [t0, t1, t2, t3, t4, t5, t6]:
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
    plt.axvline(vrv_x_layer_max, alpha=.2)
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
    plt.savefig('data/mock_vrv_av.pdf')
