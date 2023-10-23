#!/usr/bin/env python3

import argparse
import os
import re
import sys
import warnings

import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pyPLUTOplus as ppp

import cutoff_frequency


class PhaseAmplitudeData():
    def __init__(self, npz_filename):
        data = np.load(npz_filename)
        self.ini_dir = data['ini_dir'].item()
        self.vrw_x_max_layer = data['vrw_x_max_layer'].item()
        self.driver_p = data['driver_p'].item()
        self.map_t = data['map_t']
        self.map_z = data['map_z']
        self.vx1_loop_center = data['vx1_loop_center']
        self.z = data['z']
        self.dz = data['dz'].item()
        self.vp = data['vp']
        self.vp_fourier = data['vp_fourier']
        self.ck = data['ck']
        self.ωc2 = data['ωc2']
        self.ωc2_fourier = data['ωc2_fourier']

    def add_cutoff_models(self, LN_z0=u.Q(0, 'km')):
        ds = ppp.PlutoDataset(
            self.ini_dir,
            data_dir='./initial_state',
            ns_values=[0],
            )
        ds.flip_coordinate('x3', trans=ds.x3r.max())
        self.wq = cutoff_frequency.WaveQuantities(ds, args.geometry)
        self.models = [
            cutoff_frequency.CutoffSpr(self.wq),
            cutoff_frequency.CutoffBS(self.wq),
            cutoff_frequency.CutoffLN(self.wq),
            ]
        self.models_LN_offset = [
            cutoff_frequency.CutoffLN(self.wq, z0=z0) for z0 in LN_z0
            ]

    @property
    def driver_ω(self):
        return 2*np.pi/self.driver_p

    @property
    def run_id(self):
        dir_group, dir_run = self.ini_dir.strip('/').split('/')[-2:]
        reg = re.compile(r'^\d+')
        dir_group = reg.findall(dir_group)[0]
        dir_run = reg.findall(dir_run)[0]
        return f'run {dir_group}/{dir_run}'

    @property
    def labels(self):
        return {
            'period': f'$P_0 = {self.driver_p:.0f}$ s',
            'run_id': self.run_id,
            }


class PhaseAmplitudeDataCollection():
    def __init__(self, npz_filenames):
        self.data = []
        for npz_filename in npz_filenames:
            self.data.append(PhaseAmplitudeData(npz_filename))

        # if all data items share some values, set them as collection
        # attributes
        self._add_shared_quantity('vrw_x_max_layer')
        self._add_shared_quantity('map_t')
        self._add_shared_quantity('map_z')
        self._add_shared_quantity('z')
        self._add_shared_quantity('dz')
        self._add_shared_quantity('ck')

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

    def __len__(self):
        return self.data.__len__()


def compute_z_threshold(d, rv_threshold, N_hres=10000):
    rv = d.ck / d.vp
    z_hres = np.logspace(np.log10(d.z.min()), np.log10(d.z.max()), N_hres)
    rv_hres = np.interp(z_hres, d.z, rv)
    delta = np.abs(rv_hres - rv_threshold)
    mask_precision = (delta < 0.001)
    mask_vrw = (z_hres < data.vrw_x_max_layer)
    mask_derivate = (np.gradient(rv_hres) > 0)
    mask = mask_precision & mask_vrw & mask_derivate
    if np.any(mask):
        i_threshold = np.argmin(delta[mask])
        z_threshold = z_hres[mask][i_threshold]
    else:
        i_threshold = np.argmin(delta[mask_vrw])
        z_threshold = z_hres[i_threshold]
        msg = (f'threshold is below ck/vp for P0 = {d.driver_p} s '
               f'and t = {rv_threshold}. '
               f'Returning altitude of min(ck/vp).')
        warnings.warn(msg)
    return z_threshold


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('plot_dir',
                   help=('directory where to save the plots'))
    p.add_argument('phase_amplitude_data',
                   nargs='+',
                   help=('data of the runs to plot '
                         '(npz files written by phase_amplitude_3D.py)'))
    p.add_argument('-g', '--geometry', required=True,
                   help='magnetic field geometry (circular or vertical)')
    p.add_argument('--label', default='period',
                   help=('type of label to use (see options in '
                         'PhaseAmplitudeData.labels)'))
    args = p.parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)

    pdf_metadata = {'Creator': ' '.join(sys.argv)}

    data = PhaseAmplitudeDataCollection(args.phase_amplitude_data)

    for d in data:
        d.ini_dir = str(d.ini_dir)  # FIXME
        if d.ini_dir.startswith('fast_'):
            d.ini_dir = d.ini_dir[5:]
        d.add_cutoff_models(LN_z0=u.Q(np.linspace(0, 2000, 4), 'km'))

    period_colors = ['#525174', '#348AA7', '#5DD39E', '#BCE784']
    model_colors = ['#004488', '#bb5566', '#ddaa33', '#56b4e9']
    dt_res_span_kw = dict(fill=False, hatch='/', color='gray')
    vrw_span_kw = dict(fill=False, hatch='x', color='gray')

    def ylin():
        plt.yscale('linear')
    def ylog():
        plt.yscale('log')
    def xlin():
        plt.xscale('linear')
    def xlog():
        plt.xscale('log')

    with PdfPages(f'{args.plot_dir}/phase_altitude_vp.pdf', metadata=pdf_metadata) as pdf:
        plt.clf()
        for d, c in zip(data, period_colors):
            plt.plot(
                d.z, 1 / d.vp,
                '-',
                color=c,
                label=f'$1/v_p$ ({d.labels[args.label]})',
                )
            plt.plot(
                d.wq.z.to('Mm'),
                1 / d.wq.ck.to('Mm s-1'),
                '--',
                color=c,
                label=f'$1/c_k$ ({d.labels[args.label]})',
                )
        plt.legend()
        plt.axvspan(data.z.min(), data.z.min() + data.dz, **dt_res_span_kw)
        plt.axvspan(data.z.max() - data.dz, data.z.max(), **dt_res_span_kw)
        plt.axvspan(data.z.max() - data.vrw_x_max_layer, data.z.max(), **vrw_span_kw)
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$1 / v$ [s Mm⁻¹]')
        def lim():
            if plt.gca().get_yscale() == 'log':
                plt.ylim(.1, 100)
            else:
                plt.ylim(-10, 30)
            if plt.gca().get_xscale() == 'log':
                plt.xlim(4e-2, 100)
            else:
                plt.xlim(-2, 100)
        xlin(); ylin(); lim(); pdf.savefig()
        xlin(); ylog(); lim(); pdf.savefig()
        xlog(); ylin(); lim(); pdf.savefig()
        xlog(); ylog(); lim(); pdf.savefig()
        # -------------------------------
        plt.clf()
        plt.plot(
            data[0].z, 1 / data[1].vp - 1 / data[0].vp,
            'k.-',
            )
        plt.xscale('log')
        plt.xlabel('Altitude [Mm]')
        plt.ylabel(f"$\Delta(1/v_p)$")
        plt.axvspan(data.z.min(), data.z.min() + data.dz, **dt_res_span_kw)
        plt.axvspan(data.z.max() - data.dz, data.z.max(), **dt_res_span_kw)
        plt.axvspan(data.z.max() - data.vrw_x_max_layer, data.z.max(), **vrw_span_kw)
        pdf.savefig()
        # -------------------------------
        plt.clf()
        plt.plot(
            data[0].z, (1 / data[1].vp - 1 / data[0].vp) / (1 / data[0].vp),
            'k.-',
            )
        plt.xscale('log')
        plt.xlabel('Altitude [Mm]')
        plt.ylabel(f"$\Delta(1/v_p) / (1/v_{{p,\\mathrm{{{data[0].labels['run_id']}}}}})$")
        plt.axvspan(data.z.min(), data.z.min() + data.dz, **dt_res_span_kw)
        plt.axvspan(data.z.max() - data.dz, data.z.max(), **dt_res_span_kw)
        plt.axvspan(data.z.max() - data.vrw_x_max_layer, data.z.max(), **vrw_span_kw)
        pdf.savefig()

    with PdfPages(f'{args.plot_dir}/phase_altitude_vp_ck_→_ωc__debug.pdf', metadata=pdf_metadata) as pdf:
        rv_threshold = 0.4
        plt.clf()
        for d, c in zip(data, period_colors):
            z_threshold = compute_z_threshold(d, rv_threshold)
            plt.plot(d.z, d.ck / d.vp, '-', color=c, label=d.labels[args.label])
            plt.axvline(z_threshold, color=c, lw=1, zorder=1)
            plt.axhline(rv_threshold, color='k', lw=1, zorder=1)
        plt.legend()
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$c_k / v_p$')
        plt.xscale('log')
        pdf.savefig()

    with PdfPages(f'{args.plot_dir}/phase_altitude_vp_ck_→_ωc.pdf', metadata=pdf_metadata) as pdf:
        rv_thresholds = [0.2, 0.3, 0.4, 0.5]
        rv_thresholds_m = ['x', '+', 'o', 'D']
        ω = []
        z = []
        for rv_threshold in rv_thresholds:
            for d in data:
                z_threshold = compute_z_threshold(d, rv_threshold)
                ω.append(d.driver_ω)
                z.append(z_threshold)
        ω = np.array(ω).reshape((len(rv_thresholds), len(data.data)))
        z = np.array(z).reshape((len(rv_thresholds), len(data.data)))

        plt.clf()
        for m, c in zip(data[0].models, model_colors):
            if m.shortname == 'LopinNagorny2017':
                for j, m in enumerate(data[0].models_LN_offset):
                    label = m.abbr + f' {m.z0.to("km"):.0f}'
                    ls = ['-', '--', '-.', ':'][j]
                    plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=c, alpha=.5, ls=ls, label=label)
            else:
                plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=c, alpha=.5, label=m.abbr)
        for m, c in zip(data[1].models, model_colors):
            if m.shortname == 'LopinNagorny2017':
                for j, m in enumerate(data[1].models_LN_offset):
                    label = m.abbr + f' {m.z0.to("km"):.0f}'
                    ls = ['-', '--', '-.', ':'][j]
                    plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=c, ls=ls, label=label)
            else:
                plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=c, label=m.abbr)
        for i, t in enumerate(rv_thresholds):
            m = rv_thresholds_m[i]
            plt.plot(
                z[i, 0],
                ω[i, 0],
                color='k',
                alpha=.5,
                linestyle='',
                marker=m,
                fillstyle='none',
                )
            plt.plot(
                z[i, 1],
                ω[i, 1],
                color='k',
                linestyle='',
                marker=m,
                fillstyle='none',
                )
        plt.axvspan(data.z.max() - data.vrw_x_max_layer, data.z.max(), **vrw_span_kw)
        handles, labels = plt.gca().get_legend_handles_labels()
        labels += [f'Sim ($t_r = {t}$)' for t in rv_thresholds]
        handles += [mpl.lines.Line2D(
                        [], [],
                        marker=m,
                        linestyle='none',
                        fillstyle='none',
                        color='k',
                        )
                   for m in rv_thresholds_m]
        plt.legend(handles, labels, loc='lower left', ncol=2, fontsize=10)
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$\omega_c$ [s⁻¹]')
        def lim():
            if plt.gca().get_yscale() == 'log':
                plt.ylim(1e-3, 1e-1)
            else:
                plt.ylim(0, .1)
            if plt.gca().get_xscale() == 'log':
                plt.xlim(5e-2, 100)
            else:
                plt.xlim(0, 100)
        xlin(); ylin(); lim(); pdf.savefig()
        xlin(); ylog(); lim(); pdf.savefig()
        xlog(); ylin(); lim(); pdf.savefig()
        xlog(); ylog(); lim(); pdf.savefig()

    if len(data) == 2:

        with PdfPages(f'{args.plot_dir}/phase_altitude_Δzc_tr.pdf', metadata=pdf_metadata) as pdf:
            rv_thresholds = np.linspace(.2, .5, 100)
            z_thresholds = np.array([
                [compute_z_threshold(d, rv) for rv in rv_thresholds]
                for d in data])
            Δ_z_threshold = z_thresholds[1] - z_thresholds[0]
            label = (f"$z_{{c,{data[1].labels['run_id'].split(' ')[1]} }} - "
                     f" z_{{c,{data[0].labels['run_id'].split(' ')[1]} }}$")
            # -------------------------------
            plt.clf()
            plt.plot(
                rv_thresholds, Δ_z_threshold,
                'k-',
                label=label,
                )
            plt.legend()
            plt.xlabel('$r_\\omega$')
            plt.ylabel('$\\Delta z_c$ [Mm]')
            pdf.savefig()
            # -------------------------------
            plt.clf()
            plt.plot(rv_thresholds, Δ_z_threshold / z_thresholds[0], 'k-')
            plt.xlabel('$r_\\omega$')
            plt.ylabel(f"$\\Delta z_c / z_{{c,{data[0].labels['run_id'].split(' ')[1]} }}$")
            pdf.savefig()

        # compute model Δzc
        def compute_model_dz(m0, m1, ωc):
            i0 = np.argmin(np.abs(m0.ωc - ωc))
            i1 = np.argmin(np.abs(m1.ωc - ωc))
            z0 = m0.z[i0]
            z1 = m1.z[i1]
            return (z1 - z0).to('km')
        assert data[0].driver_ω == data[1].driver_ω
        ωc = u.Q(data[0].driver_ω, 's-1')
        for m0, m1 in zip(data[0].models, data[1].models):
            if m0.shortname == 'LopinNagorny2017':
                for m0, m1 in zip(data[0].models_LN_offset, data[1].models_LN_offset):
                    label = m0.abbr + f' {m0.z0.to("km"):.0f}'
                    dz = compute_model_dz(m0, m1, ωc)
                    print(label, dz)
            else:
                label = m0.abbr
                dz = compute_model_dz(m0, m1, ωc)
                print(label, dz)
