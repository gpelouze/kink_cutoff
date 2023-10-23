#!/usr/bin/env python3

import contextlib
import os
import re
import warnings

import astropy.units as u
import numpy as np
import pyPLUTOplus as ppp
import scipy.signal as ssig
import scipy.optimize as sopt
import yaml


class Summaries():
    def __init__(self, ini_dir, dbl_data_dir=None, ns_values=None, last_ns=None, load_data=True):
        self.ini_dir = ini_dir
        if not os.path.isdir(self.ini_dir):
            raise NotADirectoryError(
                f'Ini dir is not a directory: {self.ini_dir}')

        self.data_dir = f'{self.ini_dir}/output_derived/'
        if not os.path.isdir(self.data_dir):
            raise NotADirectoryError(
                f'Data dir is not a directory: {self.data_dir}')


        self.ini = ppp.PlutoIni(self.ini_dir)
        self.units = ppp.PlutoUnits(self.ini_dir)
        self.definitions = ppp.PlutoDefinitions(self.ini_dir)

        # get data_dir from ini file
        self.dbl_data_dir = dbl_data_dir
        if self.dbl_data_dir is None:
            self.dbl_data_dir = self.ini['Static Grid Output'].get('output_dir', '.')
        # Handle relative data dir definitions, eg.
        # /foo/bar -> /foo/bar/
        # ./foo/bar -> {self.ini_dir}/foo/bar/
        self.dbl_data_dir = re.sub(
            r'^\.(/|$)',
            self.ini_dir + '/',
            self.dbl_data_dir,
            )
        # pyPLUTO crashes if passed w_dir option without trailing slash
        if not self.dbl_data_dir.endswith('/'):
            self.dbl_data_dir += '/'
        if not os.path.isdir(self.dbl_data_dir):
            raise NotADirectoryError(
                f'Dbl data dir is not a directory: {self.dbl_data_dir}')

        self.grid = ppp.PlutoGridReader(self.dbl_data_dir).read()

        # load nlast_info
        with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
            self.nlast_info = ppp.pp.nlast_info(w_dir=self.dbl_data_dir)

        # determine last_ns and ns_values
        if (ns_values is not None) and (last_ns is not None):
            raise ValueError('cannot set both ns_values and last_ns')
        if ns_values is not None:
            self.last_ns = None
            self.ns_values = ns_values
            # replace negative values
            sim_last_ns = self.nlast_info['nlast']
            self.ns_values = [
                sim_last_ns + ns + 1 if ns < 0
                else ns
                for ns in self.ns_values]
        else:
            if last_ns is None:
                last_ns = self.nlast_info['nlast']
            self.last_ns = last_ns
            self.ns_values = np.arange(0, self.last_ns+1)

        self.t = self.get_t(self.ini)
        self.t_cgs = self.t * self.units.time

        self._step_data= []
        if load_data:
            self.load_data()

    def get_t(self, ini):
        snapshot_dt = float(ini['Static Grid Output']['dbl'].split(' ')[0])
        return self.ns_values * snapshot_dt

    def load_data(self):
        for ns in self.ns_values:
            self._step_data.append(self.load_step_data(ns))
            self._update_time_series()

    def load_step_data(self, ns):
        filename = f'{self.data_dir}/summary.{ns:04d}.yml'
        with open(filename, 'r') as f:
            summary = yaml.safe_load(f)
        return summary

    def _update_time_series(self):
        data = {}
        ref_data = self._step_data[0]
        # This will **probably** need recursion at some point
        for k1, v in ref_data.items():
            if not isinstance(v, dict):
                data[k1] = np.array([d[k1] for d in self._step_data])
            else:
                data[k1] = {}
                for k2, v in v.items():
                    if not isinstance(v, dict):
                        data[k1][k2] = np.array([d[k1][k2] for d in self._step_data])
                    else:
                        data[k1][k2] = {}
                        for k3, v in v.items():
                            if not isinstance(v, dict):
                                data[k1][k2][k3] = np.array([d[k1][k2][k3] for d in self._step_data])
                            else:
                                data[k1][k2][k3] = {}
                                for k4, v in v.items():
                                    if not isinstance(v, dict):
                                        data[k1][k2][k3][k4] = np.array([d[k1][k2][k3][k4] for d in self._step_data])
                                    else:
                                        data[k1][k2][k3][k4] = {}
                                        for k5, v in v.items():
                                            if not isinstance(v, dict):
                                                data[k1][k2][k3][k4][k5] = np.array([d[k1][k2][k3][k4][k5] for d in self._step_data])
                                            else:
                                                raise ValueError('summary reached max nested dict depth')

        self._data = data


def LT_over_ds(T, s):
    ''' Compute the temperature scale length divided by the grid size along the
    last axis axis of `T`.
    '''
    T_mid = (T[..., 1:] + T[..., :-1]) / 2
    dT = T[..., 1:] - T[..., :-1]
    s_mid = (s[1:] + s[:-1]) / 2
    Δs = s[1:] - s[:-1]
    L_T = T_mid / np.abs(dT/Δs)
    return s_mid, L_T / Δs


def loop_center(rho, x1, threshold=0.5):
    threshold = threshold * (
        np.min(rho, axis=0) +
        np.max(rho, axis=0)
        )
    # Domain above threshold
    mask = rho > threshold
    # Gradient along x1
    # The lower edge of the loop should contain two +0.5 cells,
    # and the upper edge two -0.5 cells.
    mask_gradient = np.gradient(mask.astype(int), axis=1)
    # Check that each time step has exactly two edges
    # More than two edges can be detected (eg. if the desity dips below the
    # threshold at the center of the loop). This is still fine, because we next
    # detect the lower and upper edges of the loop.
    number_of_edges = np.sum(np.abs(mask_gradient), axis=1)
    if not np.all(number_of_edges == 2):
        warnings.warn(f'Found more than two loop edges')
    i_x1_low = np.argmax(mask_gradient, axis=1)
    i_x1_up = len(x1) - 1 - np.argmin(mask_gradient[:, ::-1], axis=1)
    x1_low = x1[i_x1_low]
    x1_up = x1[i_x1_up]
    return (x1_low + x1_up) / 2


def compute_timelag(t, a, b, apodize=True):
    a = a - np.nanmean(a)
    b = b - np.nanmean(b)
    if apodize:
        a *= np.hanning(a.size)
        b *= np.hanning(b.size)
    norm = np.sqrt(np.sum(a**2) * np.sum(b**2))
    cc = ssig.correlate(b, a)
    cc /= norm
    if np.min(t) != 0:
        t -= np.min(t)
    tmax = np.max(t)
    tcc = np.linspace(-tmax, tmax, len(cc))
    return tcc, cc


def compute_timelag_max(tcc, cc, max_lag, fit=False):
    t1 = np.argmin(np.abs(tcc + max_lag))
    t2 = np.argmin(np.abs(tcc - max_lag))
    cc = cc[t1:t2+1]
    tcc = tcc[t1:t2+1]

    if not fit:
        imax = np.argmax(cc)
        return tcc[imax], cc[imax]

    else:
        def f_func(t, A, t0, dt, B):
            return A * np.cos((t - t0) / dt) + B

        popt, pcov = sopt.curve_fit(f_func, tcc, cc, p0=[1, 0, max_lag/4, 0])
        ccmax = popt[0] + popt[3]
        tccmax = popt[1]

        tccfit = np.linspace(tcc.min(), tcc.max(), 100)
        ccfit = f_func(tccfit, *popt)

        return tccmax, ccmax, tccfit, ccfit


class Volume():
    ''' USE analysis.derived_quantities.Volume instead '''

    def __init__(self, grid, x1lim=None, x2lim=None, x3lim=None):
        self.grid = grid

        self.x1lim = self._process_lim(grid.x1, x1lim)
        self.x2lim = self._process_lim(grid.x2, x2lim)
        self.x3lim = self._process_lim(grid.x3, x3lim)

    def _process_lim(self, x, xlim):
        if xlim is None:
            xlim = (None, None)
        xmin, xmax = xlim
        if xmin is None:
            xmin = x.xL[0]
        if xmax is None:
            xmax = x.xR[-1]
        return xmin, xmax

    @property
    def x1min(self):
        return self.x1lim[0]
    @property
    def x1max(self):
        return self.x1lim[1]
    @property
    def x2min(self):
        return self.x2lim[0]
    @property
    def x2max(self):
        return self.x2lim[1]
    @property
    def x3min(self):
        return self.x3lim[0]
    @property
    def x3max(self):
        return self.x3lim[1]

    @property
    def ix1min(self):
        return np.argmin(np.abs(self.grid.x1.xL - self.x1min))
    @property
    def ix1max(self):
        return np.argmin(np.abs(self.grid.x1.xR - self.x1max))
    @property
    def ix2min(self):
        return np.argmin(np.abs(self.grid.x2.xL - self.x2min))
    @property
    def ix2max(self):
        return np.argmin(np.abs(self.grid.x2.xR - self.x2max))
    @property
    def ix3min(self):
        return np.argmin(np.abs(self.grid.x3.xL - self.x3min))
    @property
    def ix3max(self):
        return np.argmin(np.abs(self.grid.x3.xR - self.x3max))

    @property
    def V(self):
        return ((self.x1max - self.x1min) *
                (self.x2max - self.x2min) *
                (self.x3max - self.x3min))
    @property
    def N(self):
        return self.grid.x1.N * self.grid.x2.N * self.grid.x3.N
    @property
    def dV(self):
        return self.V / self.N


class Plane():
    ''' USE analysis.derived_quantities.Plane instead '''

    def __init__(self, grid, cut_axis, cut_at):
        '''
        Parameters
        ==========
        grid : pyPLUTOplus.PlutoGrid
            Simulation grid.
        cut_axis : str ('x1', 'x2', or 'x3')
            Axis normal to the plane.
        cut_at : float
            Coordinate of the plane (along cut_axis).
        '''
        self.grid = grid
        self.cut_axis = cut_axis
        self.cut_at = cut_at

        coords = {
            'x1': ('x1', 'x2', 'x3'),
            'x2': ('x2', 'x1', 'x3'),
            'x3': ('x3', 'x1', 'x2'),
            }
        self.cut_coord, self.x_coord, self.y_coord = [
            grid.__getattribute__(cn) for cn in coords[cut_axis]]

        self.cut_axis_id = int(cut_axis.strip('x')) - 1
        self.cut_at_i = np.argmin(np.abs(self.cut_coord.x - cut_at))
        self.spatial_slice = np.roll(
            (self.cut_at_i, slice(None), slice(None)),
            self.cut_axis_id)
        self.global_slice = (slice(None), ) + tuple(self.spatial_slice)

    @property
    def xmin(self):
        return self.x_coord.xL[0]
    @property
    def xmax(self):
        return self.x_coord.xR[-1]
    @property
    def ymin(self):
        return self.y_coord.xL[0]
    @property
    def ymax(self):
        return self.y_coord.xR[-1]

    @property
    def A(self):
        return ((self.xmax - self.xmin) *
                (self.ymax - self.ymin))
    @property
    def N(self):
        return self.x_coord.N * self.y_coord.N
    @property
    def dA(self):
        return self.A / self.N


class LoopEnergyFromDataset():
    ''' USE analysis.derived_quantities.DerivedQuantitiesRecipes instead '''

    def __init__(self, ds):
        self.ds = ds
        self.γ = 5/3
        self.µ = 1/2

    @property
    def rho(self):
        return self.ds['rho'] * self.ds.units.density * u.Unit('g cm-3')
    @property
    def prs(self):
        return self.ds['prs'] * self.ds.units.pressure * u.Unit('dyn cm-2')
    @property
    def T(self):
        return (self.ds['prs'] / self.ds['rho'] * self.μ) * self.ds.units.temperature * u.Unit('K')

    @property
    def vx1(self):
        return self.ds['vx1'] * self.ds.units.velocity * u.Unit('cm s-1')
    @property
    def vx2(self):
        return self.ds['vx2'] * self.ds.units.velocity * u.Unit('cm s-1')
    @property
    def vx3(self):
        return self.ds['vx3'] * self.ds.units.velocity * u.Unit('cm s-1')
    @property
    def vmag(self):
        return np.sqrt(self.vx1**2 + self.vx2**2 + self.vx3**2)

    @property
    def Bx1(self):
        return self.ds['Bx1'] * self.ds.units.magnetic_field * u.Unit('G')
    @property
    def Bx2(self):
        return self.ds['Bx2'] * self.ds.units.magnetic_field * u.Unit('G')
    @property
    def Bx3(self):
        return self.ds['Bx3'] * self.ds.units.magnetic_field * u.Unit('G')
    @property
    def Bmag(self):
        return np.sqrt(self.Bx1**2 + self.Bx2**2 * self.Bx3**2)

    @property
    def L(self):
        return self.ds.ini['Parameters']['HALF_LOOP_L'] * u.Unit('Mm')
    @property
    def gsun(self):
        return 27.4e3 * u.Q('cm s-2')
    @property
    def g(self):
        self.gsun * np.sin((np.pi*self.ds.x3) / (2*self.L))
    @property
    def Φ(self):
        2*self.gsun*self.L/np.pi * np.cos((np.pi*self.ds.x3) / (2*self.L))

    def integrate_time(self, arr, t1, t2):
        pass

    def integrate_volume(self, arr, vol):
        sub_arr = arr[:,
                      vol.ix1min:vol.ix1max+1,
                      vol.ix2min:vol.ix2max+1,
                      vol.ix3min:vol.ix3max+1]
        # assuming all cells have the same volume
        return np.sum(sub_arr, axis=(1, 2, 3)) * vol.dV

    def integrate_surface(self, vec_arr, surf):
        ''' Parameters
        vec_arr : np.ndarray of shape (nt, nx1, nx2, nx3, 3)
            Coordinates of the vector to integrate through the surface.
        surf :
            Surface through which to integrate the vector.
            For now, the only implementation of surface is Plane.
        '''
        vec_component = vec[..., surf.cut_axis_id]
        if vec_component.ndim == 3:
            vec_component_plane = vec_component[surf.spatial_slice]
        if vec_component.ndim == 4:
            vec_component_plane = vec_component[surf.global_slice]
        return (np.sum(vec_component, axis=(1, 2)) * surf.dA)

    def EK(self, vol):
        ek = self.integrate_volume(0.5 * self.rho * self.vmag**2, vol) / vol.V
        return 1/vol.V * (ek - ek[0])

    def EM(self, vol):
        em = self.integrate_volume(self.Bmag**2 / (2*c.mu0), vol) / vol.V
        return 1/vol.V * (em - em[0])

    def EI(self, vol):
        ei = self.integrate_volume(self.prs / (γ - 1), vol) / vol.V
        return 1/vol.V * (ei - ei[0])

    def EG(self, vol):
        eg = self.integrate_volume(self.rho * self.Φ, vol) / vol.V
        return 1/vol.V * (eg - eg[0])


class LoopEnergyFromSummaries():
    def __init__(self, s):
        self.s = s

        self.volumes = {
            'full_domain': Volume(s.grid),
            }

        self.i0 = 0

    def _flux_sign(self, surf_name):
        if surf_name.endswith('BEG'):
            return +1
        elif surf_name.endswith('END'):
            return -1
        else:
            return -1

    @property
    def L(self):
        return self.s.grid.x3.xR[-1] - self.s.grid.x3.xL[0]

    @property
    def t(self):
        return self.s.t - self.s.t[self.i0]
    @property
    def t_cgs(self):
        return self.s.t_cgs - self.s.t_cgs[self.i0]
    @property
    def Dt(self):
        return np.gradient(self.s.t)

    def _int_flux(self, F, vol_name, face_name):
        F_smooth = np.zeros_like(F)
        F_smooth[0] = F[0]
        F_smooth[1:] = (F[1:] + F[:-1]) / 2
        F = F_smooth
        volume = self.s._data['volumes'][vol_name]['_volume']
        area = self.s._data['volumes'][vol_name]['_faces'][face_name]['_area']
        return np.cumsum(F * self.Dt * area / volume)

    def FS(self, vol_name, face_name):
        return self.s._data['volumes'][vol_name]['_faces'][face_name]['Svec']
    def FF(self, vol_name, face_name):
        return self.s._data['volumes'][vol_name]['_faces'][face_name]['Fvec']

    def FS_int(self, vol_name, face_name):
        return self._int_flux(self.FS(vol_name, face_name), vol_name, face_name)
    def FF_int(self, vol_name, face_name):
        return self._int_flux(self.FF(vol_name, face_name), vol_name, face_name)

    def Ein(self, vol_name):
        Ein = ( - self.FS_int(vol_name, 'X3BEG')
                + self.FS_int(vol_name, 'X3END')
                - self.FF_int(vol_name, 'X3BEG')
                + self.FF_int(vol_name, 'X3END') )
        return Ein - Ein[self.i0]

    def EM(self, vol_name):
        EM = self.s._data['volumes'][vol_name]['em']
        EM = EM - ( + self.FS_int(vol_name, 'X1END')
                    - self.FS_int(vol_name, 'X1BEG')
                    + self.FS_int(vol_name, 'X2END')
                    - self.FS_int(vol_name, 'X2BEG') )
        return EM - EM[self.i0]

    def EI(self, vol_name):
        EI = self.s._data['volumes'][vol_name]['ei']
        EI = EI - ( + self.FF_int(vol_name, 'X1END')
                    - self.FF_int(vol_name, 'X1BEG')
                    + self.FF_int(vol_name, 'X2END')
                    - self.FF_int(vol_name, 'X2BEG') )
        return EI - EI[self.i0]

    def EK(self, vol_name):
        EK = self.s._data['volumes'][vol_name]['ek']
        return EK - EK[self.i0]

    def EG(self, vol_name):
        EG = self.s._data['volumes'][vol_name]['eg']
        return EG - EG[self.i0]

    def Etot(self, vol_name):
        return self.EK(vol_name) + self.EI(vol_name) + self.EM(vol_name) + self.EG(vol_name)


class CutZDataset():
    def __init__(self, ini_dir, cut_name, data_dir=None):
        ''' Time series of cuts in the domain, written by CutZAnalysis

        Parameters
        ==========
        ini_dir : str
            Directory containing the pluto.ini and definition.h files.
            (This is not necessarily the directory containing the output data,
            i.e. output_dir defined in pluto.ini.)
        cut_name : str
            Base name of the files generated by CutZAnalysis.
            (Eg: {cut_name}.list.out and {cut_name}.0000.dat)
        data_dir : str or None (default: None)
            Directory containing PLUTO out and dbl files.
            If None, it is read from pluto.ini.
        '''
        self.ini_dir = ini_dir
        self.data_dir = data_dir
        self.cut_name = cut_name

        if not os.path.isdir(self.ini_dir):
            raise NotADirectoryError(
                f'Ini dir is not a directory: {self.ini_dir}')

        self.ini = ppp.PlutoIni(self.ini_dir)
        self.units = ppp.PlutoUnits(self.ini_dir)
        self.definitions = ppp.PlutoDefinitions(self.ini_dir)

        # get data_dir from ini file
        if self.data_dir is None:
            self.data_dir = self.ini['Static Grid Output'].get('output_dir', '.')
        # Handle relative data dir definitions, eg.
        # /foo/bar -> /foo/bar/
        # ./foo/bar -> {self.ini_dir}/foo/bar/
        self.data_dir = re.sub(
            r'^\.(/|$)',
            self.ini_dir + '/',
            self.data_dir,
            )
        # pyPLUTO crashes if passed w_dir option without trailing slash
        if not self.data_dir.endswith('/'):
            self.data_dir += '/'
        if not os.path.isdir(self.data_dir):
            raise NotADirectoryError(
                f'Data dir is not a directory: {self.data_dir}')

        self.grid = ppp.PlutoGridReader(self.data_dir).read()

        # Load cut files list
        list_fname = os.path.join(self.data_dir, f'{self.cut_name}.list.out')
        list_dtype = [
            ('nfile', int),
            ('t', float),
            ('dt', float),
            ('stepNumber', int),
            ]
        self.list = np.loadtxt(list_fname, dtype=list_dtype)

        # Load cut files data
        self.varnames = self._get_var_names(self.data_dir)
        self.data = np.ndarray((len(self.list), len(self.varnames), self.grid.x3.N))
        for i, nf in enumerate(self.list['nfile']):
            this_fname = os.path.join(
                self.data_dir,
                f'{self.cut_name}.{nf:04d}.dat')
            this_data = np.fromfile(this_fname)
            this_data = this_data.reshape(len(self.varnames), self.grid.x3.N)
            self.data[i, :, :] = this_data
        self.data = self.data.swapaxes(0, 1)

        # setup plot dir
        self.plot_dir = os.path.join(self.ini_dir, 'plots')
        os.makedirs(self.plot_dir, exist_ok=True)

    def _get_var_names(self, data_dir):
        ''' Read list of variables saved by pluto from the dbl output. '''
        varfile = os.path.join(data_dir, 'dbl.out')
        with open(varfile, 'r') as f:
            line = f.readline()
        line = line.strip().split(' ')
        varnames = line[6:]
        return varnames

    def __getitem__(self, k):
        if k in self.varnames:
            return self.data[self.varnames.index(k)]
        else:
            raise KeyError(k)

    def flip_coordinate(self, trans=0):
        ''' Flip the direction of the z coordinate

        Parameters
        ==========
        trans : float (default: 0)
            Translate the z values by the given amount after flipping.

        This function flips the coordinate sign and inverts direction of
        coordinate and var arrays, and adds `trans` to the coordinate.
        (`x3` thus becomes `-x3[::-1] + trans` and eg. `rho` becomes `rho[:,
        ::-1, :]`).
        '''
        self.grid.x3.xL = - self.grid.x3.xL[::-1] + trans
        self.grid.x3.xR = - self.grid.x3.xR[::-1] + trans
        self.data = self.data[:, :, ::-1]

    def remove_spurious_points(self):
        ''' Remove spurious points (created when restarting simulation?)


        Problem description
        -------------------
        These points have the same time and stepNumber as the next one, and the
        same values as the previous one (eg. point #3 on the next graph).

          v ┤
          a ┤          5
          l ┤       4
          u ┤    2  3
          e ┤ 1
            └─┬──┬──┬──┬──┬
               StepNumber


        Cleanup algorithm
        -----------------
            - Find array indices i such that stepNumber[i] == stepNumber[i+1]
            - For all variable var, check that np.isclose(var[i], var[i-1])
            - Remove points with indice i from time list (self.list), and data
              (self.data)
        '''

        # Find spurious points indices
        ns = self.list['stepNumber']
        i_spurious = (np.where(ns[1:] == ns[:-1])[0])
        # Check that values at spurious locations as the same as previous one
        for i in i_spurious:
            is_close = np.isclose(self.data[:, i, :], self.data[:, i-1, :])
            if not np.all(is_close):
                Ns = self.list['stepNumber'][i]
                t = self.list['t'][i]
                msg = (f'Spurious cleanup: repeated point (Ns={Ns}, t={t}),'
                       f' but values do not match.')
                warnings.warn(msg)
        self.list = np.delete(self.list, i_spurious)
        self.data = np.delete(self.data, i_spurious, axis=1)

    def trim_time(self, tmin, tmax):
        ''' Trim dataset to given time bounds

        Parameters
        ==========
        tmin : float or None
            The minimum time to keep, in seconds. If None, use the minimum time
            in the dataset.
        tmax : float or None
            The maximum time to keep, in seconds. If None, use the maximum time
            in the dataset.
        '''
        if tmin is None:
            tmin = self.t_cgs.min()
        if tmax is None:
            tmax = self.t_cgs.max()
        tmin = u.Q(tmin, 's')
        imin = np.argmin(np.abs(self.t_cgs - tmin))
        tmax = u.Q(tmax, 's')
        imax = np.argmin(np.abs(self.t_cgs - tmax))
        self.list = self.list[imin:imax+1]
        self.data = self.data[:, imin:imax+1]

    @property
    def t_cgs(self):
        return self.list['t'] * u.Q(self.units.time, 's')

    @property
    def x1(self):
        return self.grid.x1.x * u.Q(self.units.length, 'cm').to('Mm')

    @property
    def x2(self):
        return self.grid.x2.x * u.Q(self.units.length, 'cm').to('Mm')

    @property
    def x3(self):
        return self.grid.x3.x * u.Q(self.units.length, 'cm').to('Mm')

    @property
    def x1r(self):
        return self.grid.x1.xr * u.Q(self.units.length, 'cm').to('Mm')

    @property
    def x2r(self):
        return self.grid.x2.xr * u.Q(self.units.length, 'cm').to('Mm')

    @property
    def x3r(self):
        return self.grid.x3.xr * u.Q(self.units.length, 'cm').to('Mm')
