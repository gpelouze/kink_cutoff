import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

fp = mpl.font_manager.FontProperties(family=['FreeSerif', 'serif'])
plt.rcParams['figure.constrained_layout.use'] = False
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'FreeSerif'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.it'] = 'serif:italic'
plt.rcParams['mathtext.bf'] = 'serif:bold'
plt.rcParams['legend.fancybox'] = False
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.formatter.limits'] = (-5, 5)
plt.rcParams['axes.formatter.useoffset'] = False

if __name__ == '__main__':

    avm = .9995
    zv = 50
    L = 100
    tmax = 47
    t13 = tmax / 3
    Nz = 1000
    Nt = 2000

    z = np.linspace(0, L, 1000)
    t = np.linspace(0, tmax, 2000)

    izv = np.argmin(np.abs(z - zv))
    it13 = np.argmin(np.abs(t - t13))

    av3D = 1 - (1 - avm) * (z - zv) / (L - zv)
    av3D[:izv+1] = 1

    av2D = 1 - (1 - av3D.reshape(-1, 1)) * (t - t13)/t13
    av2D[:, :it13+1] = 1
    av2D[:, 2*it13:] = av3D.reshape(-1, 1)

    # paper plot
    plt.figure(0, figsize=(6.4, 3), clear=True)
    n_times = 6
    times = np.linspace(t13, 2*t13, n_times)
    time_norm = plt.matplotlib.colors.Normalize(
        vmin=0,
        vmax=times.max(),
        )
    for i, t0 in enumerate(times):
        if i == 0:
            label = f'$\\alpha_v(t \\leq \\mathrm{{{t0:.1f}~ks}})$'
        elif i == n_times - 1:
            label = (f'$\\alpha_v(t \\geq \\mathrm{{{t0:.1f}~ks}})'
                     f'= \\alpha_{{v,\\mathrm{{3D}}}}$')
        else:
            label = f'$\\alpha_v(t = \\mathrm{{{t0:.1f}~ks}})$'
        plt.plot(
            z,
            av2D[:, np.argmin(np.abs(t0 - t))],
            color=plt.cm.pink_r(time_norm(t0)),
            label=label,
            )
    plt.legend(
        frameon=False,
        )
    plt.ylabel('Velocity rewrite coefficient $\\alpha_v$')
    plt.xlabel('Altitude [Mm]')
    plt.tight_layout(pad=0.3)
    plt.savefig('data/vrw_av.pdf')
