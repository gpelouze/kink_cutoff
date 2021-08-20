#!/usr/bin/env python3

import argparse
import os

import colorednoise
import numpy as np
import pyPLUTOplus as ppp
import scipy.optimize as sopt
from matplotlib import pyplot as plt


if __name__ == '__main__':
    
    p = argparse.ArgumentParser()
    p.add_argument('ini_dir')
    p.add_argument('--beta', type=float, default=1.66,
                   help=('power law exponent (default: 1.66)'))
    p.add_argument('--seed', type=int, default=1,
                   help=('numpy random seed (default: 1)'))
    args = p.parse_args()

    np.random.seed(args.seed)
    output_dir = f'{args.ini_dir}/driver_v'
    os.makedirs(output_dir, exist_ok=True)

    # Read PLUTO configuration
    pluto_ini = ppp.PlutoIni(args.ini_dir)
    pluto_units = ppp.PlutoUnitsAstropy(args.ini_dir)

    # Verify that PLUTO is configured to run with a fixed time step
    CFL_max_var = float(pluto_ini['Time']['CFL_max_var'])
    if CFL_max_var != 1:
        msg = (f'simulations must have a fixed time step'
               f' (got CFL_max_var = {CFL_max_var} instead of 1.)')
        raise ValueError(msg)

    # Determine number of timestep
    tstop = float(pluto_ini['Time']['tstop'])
    dt = float(pluto_ini['Time']['first_dt'])
    Ns_max = int(np.ceil(tstop / dt))

    # Generate velocity time series
    v = colorednoise.powerlaw_psd_gaussian(args.beta, Ns_max)

    # Generate time arrays, in number of time steps (_ns), and seconds (_sec)
    t_ns = np.arange(Ns_max)
    t_sec = t_ns * dt * pluto_units.time.to('s')
    t_t0 = t_ns * dt

    # Save data
    dat = np.stack([t_t0, v]).T
    np.savetxt(f'{output_dir}/v.txt', dat)

    # Compute PSD of generated time series
    psd = np.abs(np.fft.rfft(v) / np.sqrt(v.size))**2
    f_ns = np.fft.rfftfreq(v.size)
    f_sec = f_ns / (dt * pluto_units.time.to('s'))

    # remove 0-frequency bin from psd
    psd = psd[1:]
    f_ns = f_ns[1:]
    f_sec = f_sec[1:]

    def fit_func(f, A, b):
        return np.log(A*f**(-b))
    popt, pcov = sopt.curve_fit(fit_func, f_ns, np.log(psd))
    psd_fit = np.exp(fit_func(f_ns, *popt))
    fit_label = 'Fit PSD = {:.2g} $\\nu^{{-{:.4g}}}$'.format(*popt)
    print(fit_label)

    # Plot data
    plt.figure(0, clear=True)
    plt.plot(t_sec.to('ks'), v, label='noise')
    plt.xlabel('Time [ks]')
    plt.ylabel('Velocity normalized to variance')
    plt.savefig(f'{output_dir}/plot_v.pdf')

    plt.figure(1, clear=True)
    plt.step(f_sec.to('Hz'), psd, label='PSD', where='mid')
    plt.plot(f_sec.to('Hz'), psd_fit, label=fit_label)
    plt.legend()
    plt.loglog()
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('PSD')
    plt.savefig(f'{output_dir}/plot_psd.pdf')
