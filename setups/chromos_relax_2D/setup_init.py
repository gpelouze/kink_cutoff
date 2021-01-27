#!/usr/bin/env python

import argparse
import os
import re
import shutil
import subprocess
import sys

import jinja2
import yaml


def load_config(filename):
    with open(filename) as f:
        config = yaml.safe_load(f)
    return config


def interpolate_pluto_ini(variables):
    template = 'pluto.ini.j2'
    output = 'pluto.ini'
    with open(template) as f:
        template = jinja2.Template(
            f.read(),
            undefined=jinja2.StrictUndefined,
            )
    with open(output, 'w') as f:
        f.write(template.render(**variables))


def make_dirs(*dir_names):
    for dir_name in dir_names:
        os.makedirs(dir_name, exist_ok=True)


def update_definitions(new_definitions):
    with open('definitions.h') as f:
        lines = f.readlines()

    in_user_def_constants = False
    new_lines = []

    # user-defined constants with values in new_definitions dictionnary
    for i, line in enumerate(lines):
        if line.startswith('/* [Beg]'):
            # start of user-defined constants
            in_user_def_constants = True

        if in_user_def_constants and line.startswith('#define'):
            _, key, value = re.split(r'\s+', line.strip())
            # If a user-defined constant in definitions.h is also present in
            # new_definitions, update its value with that of new_definitions.
            if key in new_definitions:
                value = new_definitions.pop(key)
                line = f'#define  {key}  {value}\n'

        if line.startswith('/* [End]'):
            # start of user-defined constants
            in_user_def_constants = False
            # Add remaining items of new_definitions to the definitions.h
            if new_definitions:
                for key, value in new_definitions.items():
                    new_lines.append(f'#define  {key}  {value}\n')
                new_lines.append('\n')

        new_lines.append(line)  # append current line to new_lines

    # write updated definitions lines
    with open('definitions.h',  'w') as f:
        f.writelines(new_lines)


def create_dummy_makefile(arch):
    with open('makefile', 'w') as f:
        f.write(f'ARCH = {arch}\n')


def copy_setup_files(dest_dir):
    files_to_backup = [
        'pluto.ini',
        'definitions.h',
        'makefile',
        'pluto',
        ]
    for filename in files_to_backup:
        shutil.copy(filename, dest_dir)


def get_logfile(setup_config):
    return os.path.join(setup_config['pluto_ini']['log_dir'], 'pluto.log')


def apply_setup(setup_config, no_clean=False):
    interpolate_pluto_ini(setup_config['pluto_ini'])
    update_definitions(setup_config['definitions'])
    make_dirs(setup_config['pluto_ini']['output_dir'],
              setup_config['pluto_ini']['log_dir'])
    create_dummy_makefile(setup_config['arch'])
    subprocess.run(
        ['python2',
         os.path.join(os.environ['PLUTO_DIR'], 'setup.py'),
         '--auto-update',
         ],
        check=True)
    subprocess.run(['make'], check=True)
    if not no_clean:
        subprocess.run(['make', 'clean'], check=True)


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config', default='setup_config.yml',
                   help="configuration file (default: setup_config.yml)")
    p.add_argument('--no-clean', action='store_true',
                   help="don't run 'make clean' after applying a setup")
    p.add_argument('--last-step', default=None,
                   help=("last step after which to stop "
                         "(apply_hs_equil, run_hs_equil, or apply_solve)"))
    args = p.parse_args()

    last_steps = (None, 'apply_hs_equil', 'run_hs_equil', 'apply_solve')
    if args.last_step not in last_steps:
        raise ValueError(f'unsupported last step: {args.last_step}')

    config = load_config(args.config)

    apply_setup(config['init_hs_equil'], no_clean=args.no_clean)
    copy_setup_files(config['init_hs_equil']['pluto_ini']['output_dir'])
    if args.last_step == 'apply_hs_equil':
        sys.exit(0)

    with open(get_logfile(config['init_hs_equil']), 'w') as f:
        subprocess.run(
            ['./pluto', '-maxsteps', '0'],
            stdout=f,
            stderr=f,
            check=True)
    if args.last_step == 'run_hs_equil':
        sys.exit(0)

    apply_setup(config['solve'], no_clean=args.no_clean)
    if args.last_step == 'apply_solve':
        sys.exit(0)
