#!/usr/bin/env python3

import os
import subprocess
import sys

import astropy.units as u
import numpy as np
import papy.plot
import tqdm

import matplotlib as mpl
import matplotlib.animation


u.Q = u.Quantity


def add_derived_quantities(ds):
    ds.add_var('T', ds['prs'] / ds['rho'] * 0.5)


def test_ffmpeg():
    ''' Raises RuntimeError if the ffmpeg writer is not available '''
    Writer = mpl.animation.writers['ffmpeg']


def vec_streamlines(ax, vec, x1, x2, **streamplot_kw):
    grid = np.meshgrid(x1, x2)
    kw = dict(
        linewidth=.1,
        color='w',
        arrowstyle='->',
        arrowsize=0.5,
        )
    kw.update(streamplot_kw)
    return ax.streamplot(*grid, *vec, **kw)


def apply_auto_norm(norm, arr):
    pmin, pmax = 1, 99
    if norm is 'auto_Normalize':
        vmin, vmax = np.percentile(arr, [pmin, pmax])
        return mpl.colors.Normalize(vmin, vmax)
    elif norm is 'auto_SymmetricNormalize':
        vlim = np.percentile(np.abs(arr), pmax)
        return papy.plot.SymmetricNormalize(-vlim, +vlim)
    elif norm is 'auto_LogNorm':
        vmin, vmax = np.nanpercentile(np.log10(arr), [pmin, pmax])
        return mpl.colors.LogNorm(10**vmin, 10**vmax)
    else:
        return norm


def mark_vrw_layer(ax, ds, axis='x', altitude=True,
                   color='k', alpha=0.2, linewidth=0, **kwargs):
    ''' Mark location of the velocity rewrite layer on a plot

    Parameters
    ==========
    ax : matplotlib axes
        Matplotlib axes on which to plot the mark.
    ds : pyPLUTOplus.PlutoDataset or PhaseCutZDataset
        Pluto dataset containing at least:
            - .ini : pyPLUTOplus.PlutoIni,
            - .units : pyPLUTOplus.PlutoUnits,
            - .x1, .x2, .x3 : 1D coordinate arrays.
    axis : 'x' or 'y' (default: 'x')
        Axis of the plot corresponding to the position along the loop.
    altitude : bool (default: True)
        If True, assume that in the plot, the position along the loop z is the
        altitude (z = 0 at the footpoint). If False, assume that the z = 0 at
        the loop apex.
    color, alpha, linewidth, **kwargs :
        Style properties passed to ax.axvspan or ax.axhspan.
    '''
    vrw_x_max_layer = ds.ini['Parameters'].get('vrw_x_max_layer')

    if vrw_x_max_layer is not None:
        L0 = u.Q(ds.units.length, 'cm')
        vrw_x_max_layer = (float(vrw_x_max_layer) * L0).to('Mm')

        # determine coordinate along the loop
        if ds.x3.size > 1:
            z = ds.x3
        elif ds.x2.size > 1:
            z = ds.x2
        elif ds.x1.size > 1:
            z = ds.x1
        else:
            raise ValueError('cannot determine coordinate along the loop')

        # ensure that z is in astropy units
        if not hasattr(z, 'unit'):
            z = (z * L0).to('Mm')

        # determine vrw limits
        L = z.ptp()
        z_max = vrw_x_max_layer
        z_min = z[0]  # apex
        if altitude:
            # convert coordinate with origin at loop apex to altitude
            z_max = L - z_max
            z_min = L - z_min
            z_min, z_max = z_max, z_min

        # determine appropriate ax span function
        if axis == 'x':
            axspan = ax.axvspan
        elif axis == 'y':
            axspan = ax.axhspan
        else:
            raise ValueError('invalid axis value: {axis}')

        axspan(
            z_min.to('Mm').value,
            z_max.to('Mm').value,
            color=color,
            alpha=alpha,
            linewidth=linewidth,
            **kwargs,
            )


class MapMovie():
    def __init__(self, fig, arr, coordinates=None, t_titles=None, vec=None,
                 **kwargs):
        ''' Display or save a movie of maps

        Parameters
        ==========
        fig : matplotlib.figure.Figure
            Matplotlib figure to which the movie is played.
        arr : 3D array
            Movie data, with axes as (time, y, x)
        coordinates : None or 2-tuple
            Passed to papy.plot.plot_map
        t_titles : None or 1D array
            Titles for each time timestep
        vec : None or 4D array
            Vector field with axes as (time, 2, y, x).
            Used to plot streamlines if not None.
        **kwargs :
            Passed to fig.gca().imshow through papy.plot_map.

        Methods
        =======
        play() :
            Display the movie in a matplotlib window
        save(path) :
            Save a mp4 movie.
        '''

        self.fig = fig
        self.ax = None

        self.arr = arr
        self.coordinates = coordinates
        self.t_titles = t_titles

        self.vec = vec

        # Args passed to matplotlib imshow
        self.imshow_kwargs = kwargs

        if self.imshow_kwargs.get('norm') is not None:
            self.imshow_kwargs['norm'] = apply_auto_norm(
                self.imshow_kwargs['norm'],
                self.arr,
                )

        self.plot_initiated = False

    def init_plot(self):
        self.fig.clear()
        self.ax = self.fig.gca()
        self.im = papy.plot.plot_map(
            self.ax,
            self.arr[0],
            coordinates=self.coordinates,
            animated=True,
            **self.imshow_kwargs)
        if self.vec is not None:
            self.streamlines = vec_streamlines(
                self.ax,
                self.vec[0],
                *self.coordinates,
                )
        self.cbar = self.fig.colorbar(self.im)
        self.plot_initiated = True

    def update(self, i):
        self.im.set_data(self.arr[i])
        if self.t_titles:
            self.ax.set_title(self.t_titles[i])
        else:
            self.ax.set_title('Step {:03d}'.format(i))

        # update streamlines (delete previous, add new)
        if self.vec is not None:
            self.streamlines.lines.remove()
            for art in self.ax.get_children():
                if isinstance(art, mpl.patches.FancyArrowPatch):
                    art.remove()
            self.streamlines = vec_streamlines(
                self.ax,
                self.vec[i],
                *self.coordinates,
                )

        # vlim
        if 'norm' not in self.imshow_kwargs:
            try:
                vmin = self.imshow_kwargs['vmin']
            except KeyError:
                vmin = np.nanmin(self.arr[i])
            try:
                vmax = self.imshow_kwargs['vmax']
            except KeyError:
                vmax = np.nanmax(self.arr[i])
            self.im.set_clim(vmin, vmax)

        return self.im,

    def save_frames(self, filename):
        base, ext = os.path.splitext(filename)
        frames_dir = f'{base}-frames'

        # Create frames directory and generate filenames
        os.makedirs(frames_dir, exist_ok=True)
        frames_filenames = [
            os.path.join(frames_dir, f'frame-{i:05d}.png')
            for i in range(len(self.arr))]

        # Initialize plot
        if not self.plot_initiated:
            self.init_plot()

        # Save frames
        for i, frame_filename in tqdm.tqdm(enumerate(frames_filenames)):
            self.update(i)
            self.fig.savefig(
                frame_filename,
                metadata={'Creator': ' '.join(sys.argv) + 'with viz.plot_tools.MapMovie'},
                )

        # Convert frames to video
        if os.path.exists(filename):
            os.remove(filename)
        p = subprocess.run([
            'ffmpeg',
            '-i', os.path.join(frames_dir, 'frame-%05d.png'),
            filename,
            ],
            capture_output=True,
            check=True,
            )

    def play(self):
        ''' Play the movie in a matplotlib window.  '''
        if not self.plot_initiated:
            self.init_plot()
        self.anim = mpl.animation.FuncAnimation(
            self.fig,
            self.update,
            frames=len(self.arr),
            interval=50,
            )

    def save(self, filename, fps=15, bitrate=-1):
        ''' Save a mp4 movie to `filename`. '''
        try:
            self.anim
        except AttributeError:
            self.play()
        Writer = mpl.animation.writers['ffmpeg']
        writer = Writer(fps=fps, bitrate=bitrate)
        self.anim.save(filename, writer=writer)


def rainbow_cmap():
    rgb = np.array([
        (111, 76, 155),
        (96, 89, 169),
        (84, 104, 184),
        (78, 121, 197),
        (76, 138, 198),
        (78, 150, 188),
        (83, 158, 178),
        (89, 165, 169),
        (96, 171, 158),
        (105, 177, 144),
        (119, 183, 125),
        (140, 188, 104),
        (166, 190, 84),
        (190, 188, 71),
        (209, 181, 64),
        (221, 171, 60),
        (228, 156, 57),
        (231, 140, 53),
        (230, 121, 49),
        (228, 98, 45),
        (224, 72, 40),
        (219, 33, 34),
        # (184, 34, 30),
        # (149, 32, 26),
        # (114, 30, 22),
        # (82, 25, 18),
        ]) / 255
    # https://github.com/volodia99/colorblind/
    # https://personal.sron.nl/~pault/data/colourschemes.pdf
    rgb_names = ('red', 'green', 'blue')
    x = np.linspace(0, 1, len(rgb))
    cdict = {k: list(zip(x, v, v))
             for k, v in zip(rgb_names, rgb.T)}
    return mpl.colors.LinearSegmentedColormap('colorblind_rainbow', cdict)


class PlotProperties():
    default_properties = dict(
        title='',
        label='',
        unit=1,
        cmap='gray',
        norm_class=mpl.colors.Normalize,
        )

    def __init__(self, ds):
        self.param_properties = {
            'T': dict(
                title='Temperature',
                label='Temperature [MK]',
                unit=ds.units.temperature * u.Unit('K').to('MK'),
                norm_class=mpl.colors.Normalize,
                cmap=rainbow_cmap(),
                ),

            'rho': dict(
                title='Density',
                label='Density [kg m⁻³]',
                unit=ds.units.density * u.Unit('g cm-3').to('kg m-3'),
                norm_class=mpl.colors.LogNorm,
                cmap='viridis',
                ),

            'prs': dict(
                title='Pressure',
                label='Pressure [Pa]',
                unit=ds.units.pressure * u.Unit('dyn cm-2').to('Pa'),
                cmap='magma',
                ),

            'vx1': dict(
                label='$v_x$ [km s⁻¹]',
                unit=ds.units.velocity * u.Unit('cm s-1').to('km s-1'),
                norm_class=papy.plot.SymmetricNormalize,
                cmap='PRGn',
                ),
            'vx2': dict(
                label='$v_y$ [km s⁻¹]',
                unit=ds.units.velocity * u.Unit('cm s-1').to('km s-1'),
                norm_class=papy.plot.SymmetricNormalize,
                cmap='PRGn',
                ),
            'vx3': dict(
                label='$v_z$ [km s⁻¹]',
                unit=ds.units.velocity * u.Unit('cm s-1').to('km s-1'),
                norm_class=papy.plot.SymmetricNormalize,
                cmap='PRGn',
                ),

            'Bx1': dict(
                label='$B_x$ [G]',
                unit=ds.units.magnetic_field,
                norm_class=papy.plot.SymmetricNormalize,
                cmap='PRGn',
                ),
            'Bx2': dict(
                label='$B_y$ [G]',
                unit=ds.units.magnetic_field,
                norm_class=papy.plot.SymmetricNormalize,
                cmap='PRGn',
                ),
            'Bx3': dict(
                label='$B_z$ [G]',
                unit=ds.units.magnetic_field,
                ),

            'divB': dict(
                label='divB',
                ),

            }

    def __call__(self, param, prop):
        return self.param_properties[param].get(prop, self.default_properties[prop])
