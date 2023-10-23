#!/usr/bin/env python3

import argparse
import os
import sys

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pyPLUTOplus as ppp
import scipy.optimize as sopt
import etframes

import cutoff_frequency


class AmplAltData():
    def __init__(self, npz_filename):
        data = np.load(npz_filename)
        self.ini_dir = data['ini_dir'].item()
        self.driver_p = data['driver_p'].item()
        self.vrw_x_max_layer = data['vrw_x_max_layer'].item()
        self.t = data['t']
        self.z = data['z']
        self.vx1 = data['vx1']
        self.vx1_fit = data['vx1_fit']
        self.ampl_max = data['ampl_max']
        self.ampl_fit = data['ampl_fit']
        self.P_fit = data['P_fit']

    @property
    def driver_ω(self):
        return 2*np.pi/self.driver_p


class AmplAltDataCollection():
    def __init__(self, npz_filenames):
        self.data = []
        for npz_filename in npz_filenames:
            self.data.append(AmplAltData(npz_filename))

        # if all data items share some values, set them as collection
        # attributes
        self._add_shared_quantity('t')
        self._add_shared_quantity('z')
        self._add_shared_quantity('vrw_x_max_layer')

    def _add_shared_quantity(self, q):
        ''' Set quantity as collection attribute if all collection items share
        the same value for this quantity. If not, set it to None.
        '''
        self.__setattr__(q, self.data[0].__getattribute__(q))
        for d in self.data:
            if np.all(d.__getattribute__(q) != self.__getattribute__(q)):
                self.__setattr__(q, None)

    def __getitem__(self, k):
        return self.data[k]


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('plot_dir',
                   help=('directory where to save the plots'))
    p.add_argument('amplitude_altitude_cut_data',
                   nargs='+',
                   help=('data of the runs to plot '
                         '(npz files written by amplitude_altitude_cut_3D.py)'))
    p.add_argument('--ref-run',
                   help='init directory of reference run for cutoff models')
    p.add_argument('--plot-max', action='store_true',
                   help='plot amplitude determined with max |v|')
    p.add_argument('--fit', action='store_true',
                   help='fit A(z)')
    p.add_argument('-g', '--geometry', required=True,
                   help='magnetic field geometry (circular or vertical)')
    args = p.parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)

    pdf_metadata = {'Creator': ' '.join(sys.argv)}

    data = AmplAltDataCollection(args.amplitude_altitude_cut_data)

    ds = ppp.PlutoDataset(args.ref_run, data_dir='./initial_state', ns_values=[0])
    wq = cutoff_frequency.WaveQuantities(ds, args.geometry)
    models = [
        cutoff_frequency.CutoffSpr(wq),
        cutoff_frequency.CutoffPetr(wq),
        cutoff_frequency.CutoffLN(wq),
        cutoff_frequency.CutoffBS(wq),
        ]

    period_colors = ['#525174', '#348AA7', '#5DD39E', '#BCE784']
    model_colors = ['#004488', '#bb5566', '#ddaa33', '#56b4e9']
    vrw_span_kw = dict(fill=False, hatch='x', color='gray')

    def ylin():
        plt.yscale('linear')
    def ylog():
        plt.yscale('log')
    def xlin():
        plt.xscale('linear')
    def xlog():
        plt.xscale('log')

    with PdfPages(f'{args.plot_dir}/ampl_alt_cut_vx1_ampl.pdf', metadata=pdf_metadata) as pdf:
        plt.clf()
        for d, c in zip(data, period_colors):
            plt.plot(d.z, d.ampl_fit[::-1], '-', color=c, label=f'$P_0 = {d.driver_p:.0f}$ s')
            if args.fit:
                def fit_func(z, A):
                    return A*(z-np.log(d.z.min()/2)) + np.log(2)
                popt, pcov = sopt.curve_fit(
                    fit_func,
                    np.log(d.z),
                    np.log(d.ampl_fit[::-1]),
                    )
                label = f'fit $a = {popt[0]:.3g}$'
                plt.plot(d.z, np.exp(fit_func(np.log(d.z), *popt)), '--', color=c, label=label)
            if args.plot_max:
                plt.plot(d.z, d.ampl_max[::-1], '--', color=c)
        plt.axvspan(data.z.max() - data.vrw_x_max_layer, data.z.max(), **vrw_span_kw)
        plt.legend()
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('Velocity amplitude [km s⁻¹]')
        def lim():
            if plt.gca().get_yscale() == 'log':
                plt.ylim(1, 25)
            else:
                plt.ylim(1, 22)
            if plt.gca().get_xscale() == 'log':
                plt.xlim(4e-2, 100)
            else:
                plt.xlim(-2, 100)
        xlin(); ylin(); lim(); pdf.savefig()
        xlin(); ylog(); lim(); pdf.savefig()
        xlog(); ylin(); lim(); pdf.savefig()
        xlog(); ylog(); lim(); pdf.savefig()

    with PdfPages(f'{args.plot_dir}/ampl_alt_cut_ωc.pdf', metadata=pdf_metadata) as pdf:
        plt.clf()
        for d, c in zip(data, period_colors):
            ω = u.Q(d.driver_ω, 's-1')
            z = u.Q(d.z, 'Mm')
            ck = np.interp(d.z, wq.z[::-1].to('Mm').value, wq.ck[::-1].value) * wq.ck.unit
            f = np.interp(d.z, wq.z[::-1].to('Mm').value, wq.f[::-1]) * wq.f.unit
            G = u.Q(d.ampl_fit[::-1], 'km s-1') / ω * np.sqrt(f)
            dGdz = np.gradient(G.value, z.value) * G.unit / z.unit
            d2Gdz2 = np.gradient(dGdz.value, z.value) * dGdz.unit / z.unit
            ωcSq = ω**2 + ck**2 * d2Gdz2/G
            # ωcSq_4z= ω**2 + ck**2 * (d2Gdz2/G - 1/(4*z**2))
            plt.plot(d.z, ωcSq.to('s-2'), '-', lw=2, color=c, label=f'$P_0 = {d.driver_p:.0f}$ s')
            # plt.plot(d.z, ωcSq_4z.to('s-2'), '--', lw=2, color=c)
            plt.axhline(ω.to('s-1').value**2, linestyle='--', lw=1, color=c)
        for m, c in zip(models, model_colors):
            plt.plot(m.z.to('Mm'), m.ωc**2, '-', lw=1.5, color=c, label=m.abbr)
        plt.axvspan(data.z.max() - data.vrw_x_max_layer, data.z.max(), **vrw_span_kw)
        plt.axhline(0, linewidth=1, zorder=1)
        plt.legend(ncol=2, loc='upper center')
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$\omega_c^2$ [s⁻²]')
        def lim():
            if plt.gca().get_yscale() == 'log':
                plt.ylim(1e-6, 1e-2)
            else:
                plt.ylim(-.2e-3, 1.02e-3)
            if plt.gca().get_xscale() == 'log':
                plt.xlim(5e-2, 100)
            else:
                plt.xlim(0, 100)
        xlin(); ylin(); lim(); pdf.savefig()
        xlin(); ylog(); lim(); pdf.savefig()
        xlog(); ylin(); lim(); pdf.savefig()
        xlog(); ylog(); lim(); pdf.savefig()
