#!/usr/bin/env python3


''' Plot the amplitude oscillation as a function of the position along the
loop.
'''


import argparse
import sys

import astropy.units as u
import numpy as np
import papy.plot

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pyPLUTOplus as ppp

import analysis_tools as at
import plot_tools as pt
import cutoff_frequency


# reload local packages
from importlib import reload
reload(at)
reload(pt)
reload(cutoff_frequency)


class PhaseCutZDataset(at.CutZDataset):
    def __init__(self, *args, **kwargs):
        ''' Extend CutZDataset with variables needed to compute the phase
        difference. '''
        super().__init__(*args, **kwargs)

        # Exctract data
        self.driver_p = self.ini['Parameters'].get('driver_p')
        if self.driver_p is not None:
            self.driver_p = float(self.driver_p)
            self.driver_ω = 2*np.pi/u.Q(self.driver_p, 's')
        self.vx1_loop_center = self['vx1'] * u.Q(self.units.velocity, 'cm s-1').to('km s-1').value
        self.t = self.t_cgs.to('s').value
        self.z = self.x3.to('Mm').value

        # Flip z axis
        self.vx1_loop_center = self.vx1_loop_center[:, ::-1]

class TimeShift():
    def __init__(self, dc, z0, dz):
        ''' Compute time shift between different altitudes

        Parameters
        ==========
        dc : PhaseCutZDataset
            Dataset containing velocity as a function of altitude and time.
        z0 : float
            Altitude at which to compute the time shift.
        dz : float
            Width of the window over which the time shift is computed.
        fit_cc : bool (default: True)
            If True, fit cross-correlation to find the position of its
            maximum. A cosine function is used.

        The time shift (.tshift) and phase should (.vp) should be computed and
        set by children class.
        '''
        self.z0 = z0
        self.dz = dz
        iz0 = np.argmin(np.abs(dc.z - self.z0))  # center
        izb = np.argmin(np.abs(dc.z - self.z0 + self.dz/2))  # before
        iza = np.argmin(np.abs(dc.z - self.z0 - self.dz/2))  # after
        if iza == izb:
            raise ValueError('Δz lower than grid step')
        self.vb = dc.vx1_loop_center[:, izb]
        self.va = dc.vx1_loop_center[:, iza]
        self.zb = dc.z[izb]
        self.za = dc.z[iza]

class CCTimeShift(TimeShift):
    def __init__(self, dc, z0, dz, fit_cc=False):
        ''' Compute time shift using cross-correlation

        Parameters
        ==========
        dc : PhaseCutZDataset
            Dataset containing velocity as a function of altitude and time.
        z0 : float
            Altitude at which to compute the time shift.
        dz : float
            Width of the window over which the time shift is computed.
        fit_cc : bool (default: True)
            If True, fit cross-correlation to find the position of its
            maximum. A cosine function is used.

        The time shift at altitude z0 is determined by
        computing the cross-correlation of vx1(z-dz/2) and vx1(z+dz/2) for
        different time shifts, and finding the time shift that maximizes the
        cross-correlation.
        It is then stored in `self.tshift`

        Notes:
        - This is similar to performing a time-lags analysis between
          consecutive altitude slices.
        - We restrict the search for the cross-correlation maximum to the
          interval [-P0/2, +P0/2]. In practice, the absolute cross-correlation
          maximum falls into this interval.

        '''
        super().__init__(dc, z0, dz)

        self.tcc, self.cc = at.compute_timelag(dc.t, self.vb, self.va)
        if dc.driver_p:
            self.max_lag = dc.driver_p / 4
        else:
            self.max_lag = 100  # s
        if fit_cc:
            self.tshift, self.ccmax, self.fit_tcc, self.fit_cc = at.compute_timelag_max(
                self.tcc,
                self.cc,
                self.max_lag,
                fit=True,
                )
        else:
            self.tshift, self.ccmax = at.compute_timelag_max(
                self.tcc,
                self.cc,
                self.max_lag,
                fit=False,
                )
            self.fit_tcc = None
            self.fit_cc = None

        self.vp = (self.za - self.zb) / self.tshift


class FourierTimeShift(TimeShift):
    def __init__(self, dc, z0, dz):
        ''' Compute time shift using Fourier transform

        Parameters
        ==========
        dc : PhaseCutZDataset
            Dataset containing velocity as a function of altitude and time.
        z0 : float
            Altitude at which to compute the time shift.
        dz : float
            Width of the window over which the time shift is computed.

        The time shift at altitude z0 is determined by computing the fourier
        phase at z0-dz/2 at z0+dz/2.
        '''
        super().__init__(dc, z0, dz)

        # FFT frequencies
        dt = np.gradient(dc.t)
        dt = papy.num.almost_identical(dt, .01)
        self.freqs = np.fft.rfftfreq(dc.t.size, dt)
        # FFT phase
        self.Aa, self.Pha = self._fft(self.va, dc)
        self.Ab, self.Phb = self._fft(self.vb, dc)

        imax = np.argmax(self.Aa)
        self.freq_max = self.freqs[imax]
        Pha_max = self.Pha[imax]
        Phb_max = self.Phb[imax]
        dPh = Pha_max - Phb_max
        self.tshift = - dPh / dc.driver_ω.to('s-1').value
        self.vp = (self.za - self.zb) / self.tshift

    def _fft(self, v, dc):
        F = np.fft.rfft(v * np.hanning(v.size))
        return np.abs(F), np.arctan2(F.imag, F.real)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('ini_dir',
                   help='directory containing pluto.ini')
    p.add_argument('--data-dir', default='./output',
                   help='value of output_dir from pluto.ini')
    p.add_argument('--cut-name', default='cut_center',
                   help='name of filenames generated by CutZAnalysis')
    p.add_argument('--no-reload', action='store_true',
                   help="don't load the data (for use with %run -i in IPython)")
    p.add_argument('--fit-cc', action='store_true',
                   help='fit the cross-correlation to find its maximum')
    p.add_argument('--fourier', action='store_true',
                   help='measure the phase in the Fourier domain')
    p.add_argument('--compute-cutoff-frequency', action='store_true',
                   help='derive the cutoff frequency from the phase speed')
    p.add_argument('--z0-eg', type=float, default=None,
                   help='plot phase determination example for the given z0')
    p.add_argument('--dz', type=float, default=5,
                   help='altitude step used to compute phase difference')
    p.add_argument('-g', '--geometry', required=True,
                   help='magnetic field geometry (circular or vertical)')
    args = p.parse_args()

    pdf_metadata = {'Creator': ' '.join(sys.argv)}

    if not args.no_reload:
        dc = PhaseCutZDataset(
            args.ini_dir,
            args.cut_name,
            data_dir=args.data_dir,
            )
        ds = ppp.PlutoDataset(
            args.ini_dir,
            data_dir='./initial_state',
            ns_values=[0],
            )
        ds.flip_coordinate('x3', trans=ds.x3r.max())

    if not dc.driver_p:
        if args.fourier:
            msg = 'cannot use Fourier transform with multiperiodic driver'
            raise ValueError(msg)
        if args.compute_cutoff_frequency:
            msg = 'cannot compute cutoff frequency with multiperiodic driver'
            raise ValueError(msg)

    # Cut data in time
    if dc.driver_p:
        if dc.driver_p >= 2000:
            n_periods = 1
            t_end = 2900  # s
        else:
            n_periods = 3
            t_end = dc.t.max()
        t_beg = t_end - n_periods * dc.driver_p
        i_t_beg = np.argmin(np.abs(dc.t - t_beg))
        i_t_end = np.argmin(np.abs(dc.t - t_end))
        dc.vx1_loop_center = dc.vx1_loop_center[i_t_beg:i_t_end+1]
        dc.t = dc.t[i_t_beg:i_t_end+1]

    dc.vrw_x_max_layer = float(dc.ini['Parameters'].get('vrw_x_max_layer'))

    # remove duplicate time entries
    dt = dc.t[1:] - dc.t[:-1]
    dt = np.append(dt, dt.mean())
    dt_threshold = np.mean(dt) / 2
    m = (dt > dt_threshold)
    dc.vx1_loop_center = dc.vx1_loop_center[m]
    dc.t = dc.t[m]

    # # FIXME: test with plane wave
    # dc.vx1_loop_center = 2 * np.cos(
        # dc.driver_ω.to('s-1').value * (dc.t- dc.z.reshape(-1, 1) / 1)).T

    # Compute phase difference
    z = np.linspace(dc.z.min(), dc.z.max(), int(dc.z.ptp() / args.dz))
    # Cross-correlation
    vp = np.array([
        CCTimeShift(dc, zz, args.dz, fit_cc=args.fit_cc).vp
        for zz in z
        ])
    if args.z0_eg:
        tshift_eg = CCTimeShift(dc, args.z0_eg, args.dz, fit_cc=args.fit_cc)
    # Fourier transform
    if args.fourier:
        vp_fourier = np.array([
            FourierTimeShift(dc, zz, args.dz).vp
            for zz in z
            ])
        freqs_max = np.array([
            FourierTimeShift(dc, zz, args.dz).freq_max
            for zz in z
            ])
        print('Max power periods % z:',
              'mean:', (1 / freqs_max).mean(),
              'std:', (1 / freqs_max).std(),
              )
        if args.z0_eg:
            tshift_fourier_eg = FourierTimeShift(dc, args.z0_eg, args.dz)

    # Get wave quantities
    wq = cutoff_frequency.WaveQuantities(ds, args.geometry)
    # Interpolate kink speed onto the simulation grid
    ck = np.interp(z, wq.z.to('Mm').value, wq.ck.to('Mm s-1').value)
    # Compute cutoff models
    models = [
        cutoff_frequency.CutoffSpr(wq),
        cutoff_frequency.CutoffPetr(wq),
        cutoff_frequency.CutoffLN(wq),
        cutoff_frequency.CutoffBS(wq),
        ]

    # Compute cutoff frequency from simulation data
    if args.compute_cutoff_frequency:
        ωc2 = (1 - ck**2 / vp**2) * dc.driver_ω**2
        if args.fourier:
            ωc2_fourier = (1 - ck**2 / vp_fourier**2) * dc.driver_ω**2

    # Save data for common plot
    data = dict(
        ini_dir=dc.ini_dir,
        vrw_x_max_layer=dc.vrw_x_max_layer,
        map_t=dc.t,
        map_z=dc.z,
        vx1_loop_center=dc.vx1_loop_center,
        z=z,
        dz=args.dz,
        vp=vp,
        ck=ck,
        )
    if dc.driver_p:
        data['driver_p'] = dc.driver_p
    if args.fourier:
        data['vp_fourier'] = vp_fourier
    if args.compute_cutoff_frequency:
        data['ωc2'] = ωc2
    if args.compute_cutoff_frequency and args.fourier:
        data['ωc2_fourier'] = ωc2_fourier
    np.savez(f'{dc.plot_dir}/phase_amplitude_data.npz', **data)

    plt.figure(1, clear=True)
    im = papy.plot.plot_map(
        plt.gca(),
        dc.vx1_loop_center.T,
        coordinates=(dc.t, dc.z),
        regularity_threshold=1.05,
        aspect=(dc.t.ptp() / dc.z.ptp()),
        cmap='PRGn',
        norm=papy.plot.SymmetricNormalize(),
        # vmin=-20, vmax=+20,
        )
    cb = plt.colorbar(mappable=im)
    pt.mark_vrw_layer(plt.gca(), dc, axis='y',
                      fill=False, hatch='x', alpha=.3, linewidth=1)
    if args.z0_eg:
        plt.axhline(tshift_eg.zb, color='C0')
        plt.axhline(tshift_eg.za, color='C1')
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [Mm]')
    cb.set_label('Transverse velocity [km s ⁻¹]')
    plt.savefig(f'{dc.plot_dir}/phase_altitude_vx1_t_x3.pdf', metadata=pdf_metadata)

    if args.z0_eg:
        plt.clf()
        plt.plot(dc.t, tshift_eg.vb, label=f'{tshift_eg.zb:.1f} Mm')
        plt.plot(dc.t, tshift_eg.va, label=f'{tshift_eg.za:.1f} Mm')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Transverse velocity [km s ⁻¹]')
        plt.savefig(f'{dc.plot_dir}/phase_altitude_eg_time-series.pdf', metadata=pdf_metadata)

    if args.z0_eg:
        plt.clf()
        label = f'{tshift_eg.za:.1f}–{tshift_eg.zb:.1f} Mm'
        plt.plot(tshift_eg.tcc, tshift_eg.cc, 'k-', label=label)
        plt.plot([tshift_eg.tshift], [tshift_eg.ccmax], 'ko')
        if args.fit_cc:
            plt.plot(tshift_eg.fit_tcc, tshift_eg.fit_cc, 'gray')
            plt.axvline(-tshift_eg.max_lag, linewidth=1, color='k')
            plt.axvline(+tshift_eg.max_lag, linewidth=1, color='k')
        plt.xlim(-tshift_eg.max_lag, +tshift_eg.max_lag)
        plt.axvline(0, lw=1, ls='--')
        plt.legend()
        plt.xlabel('Time delay [s]')
        plt.ylabel('Cross-correlation')
        plt.savefig(f'{dc.plot_dir}/phase_altitude_eg_cross-corr.pdf', metadata=pdf_metadata)

    if args.z0_eg and args.fourier:
        imax = np.argmax(tshift_fourier_eg.Aa)
        plt.clf()
        label = f'{tshift_eg.za:.1f}–{tshift_eg.zb:.1f} Mm'
        plt.clf()
        plt.step(tshift_fourier_eg.freqs * 1e3, tshift_fourier_eg.Aa, 'C0', where='mid', )
        # plt.step(tshift_fourier_eg.freqs * 1e3, tshift_fourier_eg.Ab, 'C1', where='mid', )
        plt.axvline(tshift_fourier_eg.freqs[imax] * 1e3, color='gray')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Frequency [mHz]')
        plt.ylabel('Normalized PSD')
        plt.savefig(f'{dc.plot_dir}/phase_altitude_eg_fourier_psd.pdf', metadata=pdf_metadata)
        plt.clf()
        plt.step(tshift_fourier_eg.freqs * 1e3, tshift_fourier_eg.Pha, 'C0', where='mid', )
        plt.step(tshift_fourier_eg.freqs * 1e3, tshift_fourier_eg.Phb, 'C1', where='mid', )
        plt.axvline(tshift_fourier_eg.freqs[imax] * 1e3, color='gray')
        plt.xscale('log')
        plt.xlabel('Frequency [mHz]')
        plt.ylabel('Phase [rad]')
        plt.savefig(f'{dc.plot_dir}/phase_altitude_eg_fourier_phase.pdf', metadata=pdf_metadata)

else:
    fft_color = '#e67931'

    def ylim_sim():
        xmax = dc.vrw_x_max_layer
        y_min = []
        y_max = []
        for line in plt.gca().lines[:1]:
            x, y = line.get_data()
            if hasattr(x, 'value'):
                x = x.value
            x = np.array(x)
            y = np.array(y)
            y = y[(y > 0) & (x < xmax)]
            if y.size > 0:
                y_min.append(np.nanmin(y))
                y_max.append(np.nanmax(y))
        if len(y_min) > 0:
            y_min = np.nanmin(y_min)
        else:
            y_min, _ = plt.ylim()
        if len(y_max) > 0:
            y_max = np.nanmax(y_max)
        else:
            _, y_max = plt.ylim()
        if plt.gca().get_yscale() == 'log':
            dy = np.log10(y_max) - np.log10(y_min)
            d = 10**(dy * .1)
            y_min /= d
            y_max *= d
        if plt.gca().get_yscale() == 'linear':
            dy = y_max - y_min
            d = dy * .1
            y_min -= d
            y_max += d
        plt.ylim(y_min, y_max)
    def ylin():
        plt.yscale('linear')
        plt.autoscale('y')
    def ylog():
        plt.yscale('log')
        plt.autoscale('y')
    def xlin():
        plt.xscale('linear')
    def xlog():
        plt.xscale('log')

    with PdfPages(f'{dc.plot_dir}/phase_altitude_vp.pdf', metadata=pdf_metadata) as pdf:
        plt.clf()
        plt.plot(z, 1 / vp, 'k.-', label='Simulation CC')
        if args.fourier:
            plt.plot(z, 1 / vp_fourier, color=fft_color, label='Simulation FFT')
        if (not args.z0_eg) and dc.driver_p:
            for i, m in enumerate(models):
                plt.plot(
                    m.z.to('Mm'),
                    1 / m.vp(dc.driver_ω).to('Mm s-1'),
                    f'C{i}-',
                    label=m.abbr)
        plt.legend()
        span_kw = dict(fill=False, hatch='/', color='gray')
        if args.z0_eg:
            plt.axvline(tshift_eg.zb, color='C0')
            plt.axvline(tshift_eg.za, color='C1')
        plt.axvspan(z.min(), z.min() + args.dz, **span_kw)
        plt.axvspan(z.max() - args.dz, z.max(), **span_kw)
        pt.mark_vrw_layer(plt.gca(), dc, axis='x',
                          fill=False, hatch='x', alpha=.3, linewidth=1)
        plt.xlabel('Altitude [Mm]')
        plt.ylabel('$1 / v_p$ [s Mm⁻¹]')
        xlin(); ylin(); ylim_sim(); pdf.savefig()
        xlin(); ylog(); ylim_sim(); pdf.savefig()
        xlog(); ylin(); ylim_sim(); pdf.savefig()
        xlog(); ylog(); ylim_sim(); pdf.savefig()

    if args.compute_cutoff_frequency:
        with PdfPages(f'{dc.plot_dir}/phase_altitude_ωc.pdf', metadata=pdf_metadata) as pdf:
            plt.clf()
            plt.plot(z, ωc2, 'k.-', label='Simulation CC')
            if args.fourier:
                plt.plot(z, ωc2_fourier, color=fft_color, label='Simulation FFT')
            for i, m in enumerate(models):
                plt.plot(m.z.to('Mm'), m.ωc, f'C{i}-', label=m.abbr)
            plt.legend()
            span_kw = dict(fill=False, hatch='/', color='gray')
            plt.axvspan(z.min(), z.min() + args.dz, **span_kw)
            plt.axvspan(z.max() - args.dz, z.max(), **span_kw)
            pt.mark_vrw_layer(plt.gca(), dc, axis='x',
                              fill=False, hatch='x', alpha=.3, linewidth=1)
            plt.xlabel('Altitude [Mm]')
            plt.ylabel('$\omega_c^2 [s⁻¹]$')
            xlin(); ylin(); ylim_sim(); pdf.savefig()
            xlin(); ylog(); ylim_sim(); pdf.savefig()
            xlog(); ylin(); ylim_sim(); pdf.savefig()
            xlog(); ylog(); ylim_sim(); pdf.savefig()
