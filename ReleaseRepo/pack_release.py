#!/usr/bin/env python3


from __future__ import print_function

import sys
from pathlib import Path
from shutil import copyfile, copytree, make_archive, move, rmtree

import f90nml

# locate base path
path_base = Path('.').resolve()

# locate build exe
path_build = path_base / 'build'
if not path_build.exists():
    print(f'{path_build} not existing! Packing stopped!')
    sys.exit()

# load version
path_sys = list(path_build.glob('*'))[0]
name_sys = path_sys.stem
path_exe = list(path_sys.glob('SUEWS*'))[0]
name_exe = path_exe.stem
name_ver = name_exe.split('_V')[-1]

# create path for archive
path_archive = path_base / ('_'.join(['SUEWS', name_ver, name_sys]))

# copy input tables
path_input_tables = path_base / 'InputTables' / name_ver
if path_input_tables.exists():
    if path_archive.exists():
        rmtree(path_archive)
    copytree(path_input_tables, path_archive)

# load path info for input and output
dict_runcontrol = f90nml.read(path_archive / 'RunControl.nml')['runcontrol']

# make input dir
path_input = dict_runcontrol['fileinputpath']
path_input = path_archive / path_input
if not path_input.exists():
    path_input.mkdir()

# move input files
list_input_files = [
    file for file in path_archive.glob('*.*')
    if 'RunControl' not in str(file)]
for file in list_input_files:
    move(str(file), str(path_input))

# make output dir
path_output = dict_runcontrol['fileoutputpath']
path_output = path_archive / path_output
if not path_output.exists():
    path_output.mkdir()

# copy SUEWS exe
copyfile(path_exe, path_archive / path_exe.name)

# archive folder as zip
make_archive(path_archive, 'zip', path_archive)

# clean workspace
rmtree(path_archive)
