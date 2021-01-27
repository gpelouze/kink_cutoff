#!/usr/bin/env python3

import numpy as np
import pyPLUTOplus as ppp
import astropy.units as u
import astropy.constants as c

import matplotlib.pyplot as plt


u.Q = u.Quantity


if __name__ == '__main__':

    ds = ppp.PlutoDataset(
        # 'data/7_chromos_relax_2D/82_visc_ramp_up_Re1e-1__300__Re1e5__1200__half_diameter/',
        'data/7_chromos_relax_2D/87_visc_ramp_up_Re1e-1__300__Re1e5__1200__half_diameter__B0=42G/',
        ns_values=[600],
        )

    # Extract values and apply normalization
    L0 = u.Q(ds.units.length, 'cm')
    L = float(ds.ini['Parameters']['half_loop_l']) * L0
    s = ds.x2 * L0
    π = u.Q(np.pi, 'rad')
    z = 2*L/np.pi * np.cos(π*s / (2*L))

    B0 = u.Q(ds.units.magnetic_field, 'G')
    Bx = ds['Bx1'][0] * B0
    Bz = ds['Bx2'][0] * B0
    B = np.sqrt(Bx**2 + Bz**2)

    rho0 = u.Q(ds.units.density, 'g cm-3')
    ρ = ds['rho'][0] * rho0

    prs0 = u.Q(ds.units.pressure, 'dyn cm-2')
    p = ds['prs'][0] * prs0

    # Cut indices
    ii = 0
    ie = ds.x1.size - 1
    i0 = ds.x2.size - 1

    # Kink speed
    Bzi = Bz[ii]
    Bze = Bz[ie]
    ρi = ρ[ii]
    ρe = ρ[ie]
    vAi = Bzi / np.sqrt(c.mu0 * ρi)
    vAe = Bze / np.sqrt(c.mu0 * ρe)
    # ck2 = (Bzi**2 + Bze**2) / (4*np.pi * (ρi + ρe))  # [51]
    ck2 = (ρi * vAi**2 + ρe * vAe**2) / (ρi + ρe)  # [51]
    ck = np.sqrt(ck2)
    ck0 = ck[i0]

    # δB2 term
    # Note: B[ii] is the same as Bz[ii], because |Bz/B - 1| < 6e-7
    B0i = B[ii, i0]
    B0e = B[ie, i0]
    δB2 = (B0i**2 - B0e**2) / (B0i**2 + B0e**2)  # [53]

    # Scale height
    pi = p[ii]
    p0i = pi[i0]
    H = - 1 / (np.gradient(np.log(pi / p0i), z.value) / z.unit)
    H0 = H[i0]
    dHdz = np.gradient(H, z.value) / z.unit

    # Cutoff frequency and period Lopin & Nagorny 2017
    ωcLN2 = ck0**2 / (4*H0*H) * (δB2 * dHdz + H**2/z**2)  # [57]
    ωcLN = np.sqrt(ωcLN2)
    PcLN = 2*np.pi / ωcLN

    # Cutoff frequency and period Spruit 1981
    g0 = u.Q(274, 'm s-2')
    g = + g0 * np.sin(π*s / (2*L))
    ωcSpr = np.sqrt(g0 / (8*H))
    PcSpr = 2*np.pi / ωcSpr
    ωcSprGProj = np.sqrt(g / (8*H))
    PcSprGProj = 2*np.pi / ωcSprGProj

    # Cutoff frequency and period Snow NAM2017
    ωcBS = vAi / (4*z)
    PcBS = 2*np.pi / ωcBS

    # Cutoff frequency Petruchkin et al. 2015
    κ = z.max() / H
    PcPetr = (np.pi**2 * κ * H) / (ck * np.sqrt(np.exp(κ) - np.exp(κ/2)))
    ωcPetr = 2*np.pi / PcPetr

    # Loop kink period
    Pk = u.Q(335, 's')
    ωk = 2*np.pi / Pk
    kwk = dict(
        color='gray',
        lw=1,
        ls='--',
        )

    # Plot
    plt.clf()
    plt.title('Pressure scale height')
    plt.plot(z.to('Mm'), H.to('m'))
    plt.yscale('log')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$H$ [m]')
    plt.savefig('data/cfc_H.pdf')

    # Plot
    plt.clf()
    plt.title('Alfvén speed')
    plt.plot(z.to('Mm'), vAi.to('m s-1'), 'C1', label='Interior')
    plt.plot(z.to('Mm'), vAe.to('m s-1'), 'C2', label='Exterior')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$v_A$ [m s⁻¹]')
    plt.savefig('data/cfc_vA.pdf')

    # Plot
    plt.clf()
    plt.title('Kink speed')
    plt.plot(z.to('Mm'), ck.to('m s-1'))
    plt.yscale('log')
    plt.xlabel('$z$ [Mm]')
    plt.ylabel('$c_k$ [m s⁻¹]')
    plt.savefig('data/cfc_ck.pdf')

    models = [
        ('Spruit (1981)', ωcSpr, PcSpr),
        ('Petrukhin (2015)', ωcPetr, PcPetr),
        ('Lopin & Nagorny (2017)', ωcLN, PcLN),
        ('Snow (NAM 2017)', ωcBS, PcBS),
        ]

    ic = 25
    x = z[ic:].to('Mm')
    x_label = '$z$ [Mm]'

    def add_propagation_condition(ax, msg):
        h, l = ax.get_legend_handles_labels()
        h = [None] + h
        l = [msg] + l
        ax.legend(h, l)

    plt.clf()
    plt.title('Kink wave cutoff frequency (propagation if $\\omega > \\omega_c$)')
    plt.xlabel(x_label)
    plt.ylabel('$\omega_c$ [s⁻¹]')
    plt.yscale('log')
    for label, ωc, _ in models:
        plt.plot(x, ωc[ic:].to('s-1'), label=label)
    plt.axhline(ωk.to('s-1').value, label='Standing kink frequency', **kwk)
    plt.ylim(6e-4, .6)
    plt.legend(frameon=False)
    plt.savefig('data/cfc_ωc.pdf')

    plt.clf()
    plt.title('Kink wave cutoff period (propagation if $P < P_c$)')
    plt.xlabel(x_label)
    plt.ylabel('$P_c$ [s]')
    plt.yscale('log')
    for label, _, Pc in models:
        plt.plot(x, Pc[ic:].to('s'), label=label)
    plt.axhline(Pk.to('s').value, label='Standing kink period', **kwk)
    plt.ylim(1e1, 1e4)
    plt.legend(frameon=False)
    plt.savefig('data/cfc_Pc.pdf')
