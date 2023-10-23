#!/usr/bin/env python3

import argparse
import collections
import os

import numpy as np
import pyPLUTOplus as ppp
import scipy.interpolate as sint
import yaml


class Metadata(dict):
    def __init__(self, ds):
        self._metadata = {
            'directories': {
                'ini_dir': ds.ini_dir,
                'data_dir': ds.data_dir,
                },
            'snapshots': {
                'ns': args.ns,
                'nlast': ds.nlast_info['nlast'],
                },
            'ini': {
                'Grid': dict(ds.ini['Grid']),
                'Time': dict(ds.ini['Time']),
                'Parameters': dict(ds.ini['Parameters']),
                },
            }

    def as_dict(self):
        return self._metadata


def is_half_diameter(x1):
    return np.all(x1 > 0)


def get_nx1(x1):
    if is_half_diameter(x1):
        nx1 = 2*len(x1)
    else:
        nx1 = len(x1)
        if nx1 % 2 != 0:
            raise ValueError('first dimension is not divisible by 2')
    return nx1


def get_r_theta(x1):
    nx1 = get_nx1(x1)
    if is_half_diameter(x1):
        r = x1
    else:
        r = x1[nx1//2:]
    theta = np.linspace(-np.pi, np.pi, nx1//2)
    return r, theta


def get_x_y_z(x1, x2, xmax=None):
    nx1 = get_nx1(x1)
    if is_half_diameter(x1):
        x1 = np.hstack((-x1[::-1], +x1))
    if xmax is None:
        xmax = x1.max() / np.sqrt(2)
    ymax = np.sqrt(x1.max()**2 - xmax**2)
    ixmax = np.argmin(np.abs(x1 - xmax))
    iymax = np.argmin(np.abs(x1 - ymax))
    # make sure that xj[ixjmax] <= xjmax
    if x1[ixmax] > xmax:
        ixmax -= 1
    if x1[iymax] > ymax:
        iymax -= 1
    margin_x = nx1 - ixmax - 1
    margin_y = nx1 - iymax - 1
    x = x1[margin_x:-margin_x]
    y = x1[margin_y:-margin_y]
    z = x2
    return x, y, z


def cartesian_to_cylindrical(x, y, z):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta, z


def extend_2D_to_3D(arr, x1, x2, vector_component=None, xmax=None):
    ''' Convert 2D cartesian array to 3D array through cylindrical symmetry.

    Input data should be a 2D loop along x2 and centered on x1 = 0.
    x1 values can be either:
        - all positive if the 2D simulation only contains half the loop in
          diameter;
        - symmetrical, centered on 0 and with a number of cells
          divisible by 2 if the 2D simulation contains the full loop diameter.

    Parameters
    ==========
    arr : 2D array of shape (nx1, nx2)
        2D cartesian data.
    x1 : 1D array of shape (nx1, )
        Coordinates along 1st axis.
    x2 : 1D array of shape (nx2, )
        Coordinates along 2nd axis.
    vector_component : 'x1', 'x2', 'x3', or None
        Whether the array contains vector components along axis x1, x2, x3, or
        a scalar (None). If it is a vector component, it is multiplied by
        cos(theta) or sin(theta) as needed.
    xmax : float or None
        Maximum value of x for the 3D array.
        The maximum value for y (ymax) is chosen such that
        xmax^2 + ymax^2 <= max(x1)^2
        If None, generate xmax = ymax = max(x1) / sqrt(2)
    '''
    r, theta = get_r_theta(x1)
    x, y, z = get_x_y_z(x1, x2, xmax=xmax)

    points = (r, theta, z)

    if is_half_diameter(x1):
        values = arr
    else:
        values = arr[len(x1)//2:]
    values = np.repeat(values, len(theta), axis=0)
    values = values.reshape(len(r), len(theta), len(z))

    xi_cartesian = np.meshgrid(x, y, z, indexing='ij')
    xi_cylindrical = cartesian_to_cylindrical(*xi_cartesian)
    xi_cylindrical = np.moveaxis(xi_cylindrical, 0, -1)

    arr_3D = sint.interpn(points, values, xi_cylindrical)

    if vector_component is not None:
        theta = xi_cylindrical[:, :, :, 1]
        if vector_component == 'x1':
            arr_3D *= np.cos(theta)
        elif vector_component == 'x2':
            arr_3D *= np.sin(theta)
        elif vector_component == 'x3':
            pass
        else:
            raise ValueError(f'unknown vector component: {vector_component}')

    return arr_3D


def get_vector_component(varname):
    if varname.endswith('x1'):
        return 'x1'
    elif varname.endswith('x2'):
        return 'x2'
    elif varname.endswith('x3'):
        return 'x3'
    else:
        return None


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('ini_dir')
    p.add_argument('--data-dir', default='./output')
    p.add_argument('--ns', type=int, default=-1)
    p.add_argument('--xmax', type=float)
    p.add_argument('-O', '--overwrite', action='store_true')
    args = p.parse_args()

    ds = ppp.PlutoDataset(
        args.ini_dir,
        data_dir=args.data_dir,
        ns_values=[args.ns],
        )

    output_3D_dir = f'{ds.ini_dir}/output_3D'
    try:
        os.makedirs(output_3D_dir)
    except FileExistsError as e:
        if not args.overwrite:
            msg = "3D data already exist, use '-O' to overwrite"
            raise FileExistsError(msg) from e

    metadata = Metadata(ds)
    with open(os.path.join(output_3D_dir, 'metadata.yml'), 'w') as f:
        yaml.safe_dump(metadata.as_dict(), f, sort_keys=False)

    varnames = [
        ('rho', 'rho'),
        ('vx1', 'vx1'),
        ('vx1', 'vx2'),
        ('vx2', 'vx3'),
        ('Bx1', 'Bx1'),
        ('Bx1', 'Bx2'),
        ('Bx2', 'Bx3'),
        ('prs', 'prs'),
        ('psi_glm', 'psi_glm'),
        ]
    data = collections.OrderedDict()
    for varname_2D, varname_3D in varnames:
        print(f'Converting {varname_3D}')
        arr = ds[varname_2D][0]
        arr = extend_2D_to_3D(
            arr, ds.x1, ds.x2,
            vector_component=get_vector_component(varname_3D),
            xmax=args.xmax,
            )
        arr = arr.reshape(1, *arr.shape)
        data[varname_3D] = arr
    t = np.array([0.])
    x1, x2, x3 = get_x_y_z(ds.x1, ds.x2, xmax=args.xmax)

    # save PLUTO files
    dataset = ppp.DblDataset(
        data=data,
        coordinates=(t, x1, x2, x3),
        n_dimensions=3,
        geometry='cartesian',
        )
    dataset.save_dbl(output_3D_dir, file_type='multiple_files')
