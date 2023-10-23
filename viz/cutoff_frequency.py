#!/usr/bin/env python3

import argparse
import os
import sys

import numpy as np
import pyPLUTOplus as ppp
import astropy.units as u
import astropy.constants as c
import papy.plot
import pandas as pd
import scipy.integrate
import scipy.interpolate
import tqdm

import matplotlib.pyplot as plt

import plot_tools as pt

# reload local packages
from importlib import reload
reload(pt)


u.Q = u.Quantity


class WaveQuantities():
    def __init__(self, ds, geometry):
        if ds.ndim not in (2, 3):
            raise ValueError(f'unsupported dimension: {ds.ndim}')

        # Extract values and apply normalization
        L0 = u.Q(ds.units.length, 'cm')
        self.L = float(ds.ini['Parameters']['half_loop_l']) * L0
        if ds.ndim == 2:
            self.s = ds.x2 * L0
        elif ds.ndim == 3:
            self.s = ds.x3 * L0
        π = u.Q(np.pi, 'rad')

        # Altitude and projected gravity
        g0 = u.Q(274, 'm s-2')
        if geometry == 'circular':
            self.z = 2*self.L/np.pi * np.cos(π*(self.L - self.s) / (2*self.L))
            self.g = g0 * np.sin(π*(self.L - self.s) / (2*self.L))
        elif geometry == 'vertical':
            self.z = self.s
            self.g = g0
        else:
            raise ValueError(f'unknown geometry: {geometry}')

        B0 = u.Q(ds.units.magnetic_field, 'G')
        if ds.ndim == 2:
            Bz = ds['Bx2'][0] * B0
        elif ds.ndim == 3:
            Bz = ds['Bx3'][0] * B0

        rho0 = u.Q(ds.units.density, 'g cm-3')
        ρ = ds['rho'][0] * rho0

        prs0 = u.Q(ds.units.pressure, 'dyn cm-2')
        p = ds['prs'][0] * prs0

        if ds.ndim == 3:
            # remove y dimension
            iy = np.argmin(np.abs(ds.x2 - 0.))
            Bz = Bz[:, iy]
            ρ = ρ[:, iy]
            p = p[:, iy]

        assert Bz.ndim == 2

        # Cut indices
        ii = np.argmin(np.abs(ds.x1 - 0.))
        ie = np.argmin(np.abs(ds.x1 - ds.x1.max()))
        if ds.ndim == 2:
            self.i0 = ds.x2.size - 1
        elif ds.ndim == 3:
            self.i0 = ds.x3.size - 1

        # Kink speed
        self.Bzi = np.abs(Bz[ii])
        self.Bze = np.abs(Bz[ie])
        ρi = ρ[ii]
        ρe = ρ[ie]
        self.vAi = self.Bzi / np.sqrt(c.mu0 * ρi)
        self.vAe = self.Bze / np.sqrt(c.mu0 * ρe)
        # ck2 = (Bzi**2 + Bze**2) / (4*np.pi * (ρi + ρe))  # [51]
        ck2 = (ρi * self.vAi**2 + ρe * self.vAe**2) / (ρi + ρe)  # [51]
        self.ck = np.sqrt(ck2)
        # self.ck[:] = u.Q(1e5, 'm s-1')  # FIXME: constant ck

        # Scale height
        pi = p[ii]
        p0i = pi[self.i0]
        self.H = - 1 / (np.gradient(np.log(pi / p0i), self.z.value) / self.z.unit)

        self.I = np.cumsum(1/self.H * np.gradient(self.z)).to('')
        self.f = np.exp(- 0.5 * self.I)


class VALWaveQuantities():
    def __init__(self):
        # Compatibility attrs from WaveQuantities
        # self.L = None
        # self.s = None
        self.g = u.Q(274, 'm s-2')

        # Extract values
        # Load VAL-C
        df = pd.read_csv(
            os.path.join(
                os.getenv('SOL', os.environ['HOME']),
                'simu/VAL-atmospheres/data/val-c.dat',
                ),
            comment='#',
            )
        # df = df[17:]  # FIXME
        z = u.Q(df['h'], 'km')
        ρ = c.m_p * u.Q(df['ne'], 'cm-3')
        p = u.Q(df['Pt'], 'dyn cm-2')

        # Upscale
        self.z = np.linspace(z[0], z[-1], 2048)
        ρ = self.interp(self.z, z, ρ)
        p = self.interp(self.z, z, p)

        self.ρ = ρ
        self.p = p

        # Interior and exterior quantities
        B0 = u.Q(22.8, 'G')  # FIXME: assumption
        self.Bzi = np.repeat(B0, self.z.size)
        self.Bze = np.repeat(B0, self.z.size)
        ρe = ρ
        ρi = ρ * 3  # FIXME: assumption

        # Cut indices

        # Kink speed
        self.vAi = self.Bzi / np.sqrt(c.mu0 * ρi)
        self.vAe = self.Bze / np.sqrt(c.mu0 * ρe)
        ck2 = (ρi * self.vAi**2 + ρe * self.vAe**2) / (ρi + ρe)  # [51]
        self.ck = np.sqrt(ck2)

        # Scale height
        self.i0 = np.argmin(np.abs(self.z))  # doesn’t actually change H
        p0 = p[0]
        self.H = - 1 / (np.gradient(np.log(p / p0), self.z.value) / self.z.unit)

        self.I = np.cumsum(1/self.H * np.gradient(self.z)).to('')
        self.f = np.exp(- 0.5 * self.I)

    def interp(self, new_x, x, y):
        tck = scipy.interpolate.splrep(x.value[::-1], y.value[::-1])
        return scipy.interpolate.splev(new_x, tck) * y.unit
        # return np.interp(new_x.value, x.value[::-1], y.value[::-1]) * y.unit


class CutoffModel():
    def __init__(self, wq, ωc):
        self.wq = wq
        self.ωc = ωc.to('s-1')
        self.Pc = (2*np.pi / ωc).to('s')

    def lc2(self, ω):
        return self.wq.ck**2 / (ω**2 - self.ωc**2)

    def Lc(self, ω):
        ic = 26  # FIXME
        b = ω.reshape(-1, 1) > self.ωc[:-ic]
        i = np.argmax(b, axis=1)
        Lc = self.wq.z[:-ic][i]
        # propagation for all z => no tunneling
        Lc[np.all(b, axis=1)] = np.nan
        # propagation for no z => infinite tunneling
        Lc[np.all(~b, axis=1)] = np.inf
        if np.any(np.isinf(Lc)):
            Lc[np.isinf(Lc).argmax()] = self.wq.z[:-ic].max()  # plot until domain end
        return Lc

    @property
    def fourzSq(self):
        return 4*self.z**2

    def Q(self, ω):
        Q = 1 / self.lc2(ω)
        Q = Q - 1 / self.fourzSq
        return Q

    @property
    def z(self):
        return self.wq.z

    def vp2(self, ω):
        ''' Phase speed squared '''
        return self.wq.ck**2 * ω**2 / (ω**2 - self.ωc**2)

    def vp(self, ω):
        ''' Phase speed '''
        return np.sqrt(self.vp2(ω))


class CutoffStep(CutoffModel):
    label = 'Step cutoff'
    shortname = 'StepCutoff'
    abbr = 'Step'
    def __init__(self, wq):
        ωc_low = 2*np.pi / u.Q(1000, 's')
        ωc_high = 2*np.pi / u.Q(300, 's')
        ωc = np.full(wq.z.shape, ωc_low) * ωc_low.unit
        ωc[(u.Q(20, 'Mm') < wq.z) & (wq.z < u.Q(30, 'Mm'))] = ωc_high
        super().__init__(wq, ωc)


class CutoffSpr(CutoffModel):
    label = 'Spruit (1981)'
    shortname = 'Spruit1981'
    abbr = 'Sp81'
    def __init__(self, wq):
        ωc = np.sqrt(wq.g / (8*wq.H))
        # ωcGProj = np.sqrt(g / (8*H))
        super().__init__(wq, ωc)


class CutoffPetr(CutoffModel):
    label = 'Petruchkin et al. (2015)'
    shortname = 'Petruchkin+2015'
    abbr = 'P15'
    def __init__(self, wq):
        κ = wq.z.max() / wq.H
        Pc = (np.pi**2 * κ * wq.H) / (wq.ck * np.sqrt(np.exp(κ) - np.exp(κ/2)))
        ωc = 2*np.pi / Pc
        super().__init__(wq, ωc)


class CutoffLN(CutoffModel):
    label = 'Lopin & Nagorny (2017)'
    shortname = 'LopinNagorny2017'
    abbr = 'LN17'
    def __init__(self, wq, z0=u.Q(0, 'Mm')):
        i0 = np.argmin(np.abs(wq.z - z0))
        # Note: B[ii] is the same as Bz[ii], because |Bz/B - 1| < 6e-7
        B0i = wq.Bzi[i0]
        B0e = wq.Bze[i0]
        δB2 = (B0i**2 - B0e**2) / (B0i**2 + B0e**2)  # [53]
        ck0 = wq.ck[i0]
        H0 = wq.H[i0]
        dHdz = np.gradient(wq.H, wq.z.value) / wq.z.unit
        self.z0 = wq.z[i0]

        ωc2 = ck0**2 / (4*H0*wq.H) * (δB2 * dHdz + wq.H**2/wq.z**2)  # [57]
        # ωc2 = wq.ck**2 / (4*wq.H**2) * (δB2 * dHdz + wq.H**2/wq.z**2)  # [57]
        ωc = np.sqrt(ωc2)
        super().__init__(wq, ωc)


class CutoffBS(CutoffModel):
    label = 'Snow (NAM 2017)'
    shortname = 'Snow2017'
    abbr = 'Sn17'
    def __init__(self, wq):
        ωc = wq.vAi / (4*wq.z)
        super().__init__(wq, ωc)

class CutoffRA(CutoffModel):
    ''' WARNING: this is the GLOBAL cutoff frequency, but the paper also has a
    local one '''
    label = 'Routh et al. (2013) A'
    shortname = 'Routh+2013A'
    abbr = 'R13A'
    def __init__(self, wq):
        ωc = wq.ck / (4*wq.H)
        super().__init__(wq, ωc)


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('ini_dir')
    p.add_argument('data_dir')
    p.add_argument('-n', '--ns', type=int, required=True)
    p.add_argument('-g', '--geometry', required=True, help='magnetic field geometry (circular or vertical)')
    args = p.parse_args()

    pdf_metadata = {'Creator': ' '.join(sys.argv)}

    ds = ppp.PlutoDataset(
        args.ini_dir,
        data_dir=args.data_dir,
        ns_values=[args.ns],
        )
    if ds.ndim == 2:
        ds.flip_coordinate('x2', trans=ds.x2r.max())
    elif ds.ndim == 3:
        ds.flip_coordinate('x3', trans=ds.x3r.max())

    ds.plot_dir = os.path.join(ds.ini_dir, 'plots')
    os.makedirs(ds.plot_dir, exist_ok=True)

    wq = WaveQuantities(ds, args.geometry)

    models = [
        CutoffSpr(wq),
        CutoffPetr(wq),
        CutoffLN(wq),
        CutoffBS(wq),
        # CutoffRA(wq),
        ]

    # Explored frequency range
    P_arr = np.logspace(1, 4, 500) * u.Q(1, 's')
    ω_arr = 2*np.pi / P_arr
    # Select frequencies
    P0_list = u.Q([200, 335, 700, 2000], 's')
    P0_color = ['#525174', '#348AA7', '#5DD39E', '#BCE784']
    ω0_list = 2*np.pi / P0_list
    iP0 = np.argmin(np.abs(P0_list.reshape(-1, 1) - P_arr), axis=1)
    kw0 = dict(
        color='gray',
        lw=1,
        ls='--',
        )

    plt.clf()
    plt.title('Pressure scale height')
    plt.plot(wq.z.to('Mm'), wq.H.to('m'))
    plt.yscale('log')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$H$ [m]')
    pt.mark_vrw_layer(plt.gca(), ds)
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_H.pdf', metadata=pdf_metadata)

    plt.clf()
    plt.title('Alfvén speed')
    plt.plot(wq.z.to('Mm'), wq.vAi.to('m s-1'), 'C1', label='Interior')
    plt.plot(wq.z.to('Mm'), wq.vAe.to('m s-1'), 'C2', label='Exterior')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$v_A$ [m s⁻¹]')
    pt.mark_vrw_layer(plt.gca(), ds)
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_vA.pdf', metadata=pdf_metadata)

    plt.clf()
    plt.title('Kink speed')
    plt.plot(wq.z.to('Mm'), wq.ck.to('m s-1'))
    plt.yscale('log')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$c_k$ [m s⁻¹]')
    pt.mark_vrw_layer(plt.gca(), ds)
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_ck.pdf', metadata=pdf_metadata)

    plt.clf()
    plt.title('Kink wave cutoff frequency (propagation if $\\omega > \\omega_c$)')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$\omega_c$ [s⁻¹]')
    plt.yscale('log')
    for m in models:
        plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), label=m.label)
    for ω0 in ω0_list:
        plt.axhline(ω0.to('s-1').value, **kw0)
    plt.ylim(6e-4, .6)
    plt.legend(frameon=False)
    pt.mark_vrw_layer(plt.gca(), ds)
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_ωc.pdf', metadata=pdf_metadata)
    plt.xscale('log')
    plt.savefig(f'{ds.plot_dir}/cfc_ωc_log.pdf', metadata=pdf_metadata)

    models_LN_offset = [
        CutoffLN(wq, z0=z0) for z0 in u.Q(np.linspace(0, 1000, 4), 'km')
        ]
    plt.clf()
    plt.title('Kink wave cutoff frequency (propagation if $\\omega > \\omega_c$)')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$\omega_c$ [s⁻¹]')
    plt.yscale('log')
    for i, m in enumerate(models):
        if m.shortname == 'LopinNagorny2017':
            for j, m in enumerate(models_LN_offset):
                label = m.abbr + f' {m.z0.to("km"):.0f}'
                ls = ['-', '--', '-.', ':'][j]
                plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=f'C{i}', ls=ls, label=label)
        else:
            plt.plot(m.z.to('Mm'), m.ωc.to('s-1'), color=f'C{i}', label=m.abbr)
    for ω0 in ω0_list:
        plt.axhline(ω0.to('s-1').value, **kw0)
    plt.ylim(6e-4, .6)
    plt.legend(frameon=False, loc='upper right', ncol=2)
    pt.mark_vrw_layer(plt.gca(), ds)
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_ln_ωc.pdf', metadata=pdf_metadata)
    plt.xscale('log')
    plt.savefig(f'{ds.plot_dir}/cfc_ln_ωc_log.pdf', metadata=pdf_metadata)

    plt.clf()
    plt.title('Kink wave cutoff period (propagation if $P < P_c$)')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$P_c$ [s]')
    plt.yscale('log')
    for m in models:
        plt.plot(m.z.to('Mm'), m.Pc.to('s'), label=m.label)
    for P0 in P0_list:
        plt.axhline(P0.to('s').value, **kw0)
    plt.ylim(1e1, 1e4)
    plt.legend(frameon=False)
    pt.mark_vrw_layer(plt.gca(), ds)
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_Pc.pdf', metadata=pdf_metadata)
    plt.xscale('log')
    plt.savefig(f'{ds.plot_dir}/cfc_Pc_log.pdf', metadata=pdf_metadata)

    plt.clf()
    plt.title('Tunneling length')
    plt.xlabel('P [s]')
    plt.ylabel('$L_c$ [Mm]')
    plt.xscale('log')
    plt.yscale('log')
    for i, m in enumerate(models):
        Lc = m.Lc(2*np.pi/P_arr)
        plt.plot(P_arr.to('s'), Lc.to('Mm'), label=m.label)
    for P0 in P0_list:
        plt.axvline(P0.to('s').value, **kw0)
    plt.legend(frameon=False)
    # plt.ylim(1e-3, 1e2)
    pt.mark_vrw_layer(plt.gca(), ds, axis='y')
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_Lc.pdf', metadata=pdf_metadata)

    plt.clf()
    plt.title('Damping length')
    plt.xlabel('z [s]')
    plt.ylabel('$l_c$ [Mm]')
    plt.yscale('log')
    plt.xscale('log')
    for i, m in enumerate(models):
        lc2 = m.lc2(2*np.pi / u.Q(200, 's'))
        plt.plot(m.z.to('Mm'), np.sqrt(lc2).to('Mm'), f'C{i}-', label=m.label)
        plt.plot(m.z.to('Mm'), np.sqrt(-lc2).to('Mm'), f'C{i}--')
    plt.legend(frameon=False)
    plt.ylim(1e-2, 1e3)
    pt.mark_vrw_layer(plt.gca(), ds)
    plt.tight_layout()
    plt.savefig(f'{ds.plot_dir}/cfc_lc.pdf', metadata=pdf_metadata)

    # -------------------------------------------------------------------------

    # save data for numerical resolution of ODE
    rk_zmax = 90.
    rk_dz = .01
    rk_n_step = np.ceil(rk_zmax / rk_dz) + 1
    rk_z = np.linspace(0, rk_zmax, rk_n_step)
    z = wq.z.to('Mm').value
    # save f
    f = wq.f
    rk_f = np.interp(rk_z, z, f)
    np.savetxt(f'{ds.plot_dir}/f.txt', np.stack([rk_z, rk_f.to('').value]).T)
    # save q
    for m in models:
        ax_Q = plt.figure(1, clear=True).gca()
        ax_lc = plt.figure(2, clear=True).gca()
        for i, ω0 in enumerate(ω0_list):
            Q = m.Q(ω0).to('Mm-2').value
            llc = (1/m.lc2(ω0)).to('Mm-2').value
            rk_Q = np.interp(rk_z, z, Q)
            # save data
            P0 = 2*np.pi / ω0
            P0_str = f'{P0:.0f}'.replace(' ', '')
            np.savetxt(
                f'{ds.plot_dir}/Q_{m.shortname}_{P0_str}.txt',
                np.stack([rk_z, rk_Q]).T,
                )
            # plot
            ax_Q.plot(z, Q, '-', color=P0_color[i], label=f'$P=${P0:.0f}')
            ax_Q.plot(z, -Q, '--', color=P0_color[i])
            ax_lc.plot(z, llc, '-', color=P0_color[i], label=f'$P=${P0:.0f}')
            ax_lc.plot(z, -llc, '--', color=P0_color[i])
        zz = 1/m.fourzSq.to('Mm2').value
        ax_Q.plot(z, zz, 'r-', alpha=.5, label='$-1/(4z^2)$')
        for ax in (ax_Q, ax_lc):
            ax.set_xlabel('$z$ [Mm]')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend()
        ax_Q.set_ylabel('$Q(\omega)$ [Mm⁻²]')
        ax_lc.set_ylabel('$l_c^{-2}(\omega)$ [Mm⁻²]')
        ax_Q.figure.savefig(f'{ds.plot_dir}/Q_{m.shortname}_Q.pdf')
        ax_lc.figure.savefig(f'{ds.plot_dir}/Q_{m.shortname}_lc.pdf')

    # Solve wave equation
    # In this part, z=0 is at the first cell of arrays (instead of the last
    # like above)
    sol_dir = f'{ds.plot_dir}/wave_equation'
    os.makedirs(sol_dir, exist_ok=True)
    z = wq.z.to('Mm').value
    f = wq.f
    # save q
    ξ_Lc = {}
    ξ_ωz = {}
    res_P0 = {}
    # models_to_solve = models
    models_to_solve = [models[0], models[2], models[3]]
    # models_to_solve = [models[2]]
    for m in models_to_solve:
        ξ_Lc[m] = np.full_like(ω_arr.value, -1)
        ξ_ωz[m] = np.full((ω_arr.size, z.size), -1.)
        res_P0[m] = []
        for i, ω in enumerate(tqdm.tqdm(ω_arr)):
            Q = m.Q(ω).to('Mm-2').value
            Q = scipy.interpolate.interp1d(z, Q, fill_value='extrapolate')
            def wave_equation_rhs(z, F):
                return [F[1], -Q(z) * F[0]]
            y0 = 1
            if Q(0) < 0:
                y0der = +y0 * np.sqrt(-Q(0))
            else:
                y0der = 0
            res = scipy.integrate.solve_ivp(
                wave_equation_rhs,
                t_span=(z[0], z[-1]),
                y0=[y0, y0der],
                t_eval=z,
                )
            ξ = res.y[0] / np.sqrt(f)
            ξ_ωz[m][i, :] = ξ
            Lc = m.Lc(ω).to('Mm').value[0]
            i_Lc = np.argmin(np.abs(z - Lc))
            ξ_Lc[m][i] = ξ[i_Lc] / ξ[0]
            # Plot differential equation result for frequencies in ω0_list
            if i in iP0:
                res_P0[m].append(res)

    for m in models_to_solve:
        plt.clf() # -----------------------------------------------------------
        for i, (P0, res) in enumerate(zip(P0_list, res_P0[m])):
            ξ = res.y[0] / np.sqrt(f)
            Lc = m.Lc(2*np.pi/P0).to('Mm').value[0]
            plt.plot(res.t, ξ / ξ[0], color=P0_color[i], label=f'P = {P0:.0f}')
            plt.axvline(Lc, color=P0_color[i], ls='--')
        plt.legend()
        plt.xlabel('$z$ [Mm]')
        plt.ylabel('$\\xi(z) / \\xi(0)$')
        plt.savefig(f'{sol_dir}/{m.shortname}_ξ_z.pdf')
        plt.clf() # -----------------------------------------------------------
        plt.plot(P_arr.to('s'), ξ_Lc[m])
        for P0 in P0_list:
            plt.axvline(P0.to('s').value, **kw0)
        plt.xscale('log')
        plt.xlabel('$P$ [s]')
        plt.ylabel('$\\xi(L_c) / \\xi(0)$')
        plt.savefig(f'{sol_dir}/{m.shortname}_ξ_Lc.pdf')
        plt.clf() # -----------------------------------------------------------
        x = np.log10(P_arr.to('s').value)
        im = papy.plot.plot_map(
            plt.gca(),
            (ξ_ωz[m] / ξ_ωz[m][:, 0].reshape(-1, 1)),
            coordinates=[z, x],
            aspect=(z.ptp() / x.ptp()),
            )
        plt.ylabel('log$_{10}$($P$ [s])')
        plt.xlabel('$z$ [Mm]')
        plt.colorbar(im, label='$\\xi(z) / \\xi(0)$')
        plt.savefig(f'{sol_dir}/{m.shortname}_ξ_ωz.pdf')

    Δv = np.array([  # P0 [s], Δvx1(z) / Δvx1(L)
        # FIXME: coarse values measured from plot!
        (200, 16/2),
        (335, 6.5/2),
        (700, 3.2/2),
        (2000, 2.3/2),
        ])
    plt.clf() # -----------------------------------------------------------
    z0 = 40  # Mm
    iz0 = np.argmin(np.abs(z - z0))
    for m in models_to_solve:
        ξ = ξ_ωz[m][:, iz0] / ξ_ωz[m][:, 0]
        plt.plot(P_arr.to('s'), ξ, label=m.shortname)
    plt.plot(Δv.T[0], Δv.T[1], 'kx', label='3D runs')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$P$ [s]')
    plt.ylabel('$\\xi(z0) / \\xi(0)$')
    plt.savefig(f'{sol_dir}/all_ξ_ω.pdf')
