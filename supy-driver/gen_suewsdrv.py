#!/usr/bin/env python
#####################################################################
# generator of module_sf_suewsdrv.F by combining source files
# authors:
# Zhenkun Li, lizhenk@yeah.net
# Ting Sun, ting.sun@reading.ac.uk
# history:
# 13 Aug 2018, initial version
#####################################################################

import pandas as pd
import numpy as np
from glob import glob
import os
from copy import copy


def get_file_list(path_Makefile):
    """Short summary.

    Parameters
    ----------
    path_Makefile : string
        path to makefile for dependencies

    Returns
    -------
    List
        a list of dependencies

    """
    # read in makefile source code as pd.Series
    code_raw = pd.read_csv(path_Makefile, header=None, sep='\n',
                           engine='python', comment='#', squeeze=True)
    # clean source code
    code_clean = code_raw.str.replace('\t', ' ', regex=False).str.strip()

    # retrieve lines for dependencies
    modules = [
        'UTILS', 'MODULES', 'OTHERS', 'TEST',
        # 'WRF', # we dont need 'WRF' for supy
    ]
    # positions for staring lines
    pos_file_start = [code_clean.index[code_clean.str.startswith(mod)][0]
                      for mod in modules]
    # positions for ending lines
    pos_file_end = copy(pos_file_start[1:]) + [pos_file_start[-1] + 1]

    # line blocks of groups
    lines_mod = [code_clean.iloc[start:end]
                 for start, end in zip(pos_file_start, pos_file_end)]

    # organise dependencies as dicts for groups
    mod_files = [mod.str.replace('\\', '').str.split('=').sum()
                 for mod in lines_mod]
    list_mod_files = [pd.Series(mod[1:]).str.strip() for mod in mod_files]

    # combine all files into one list
    list_files = pd.concat(list_mod_files).reset_index(
        drop=True).str.replace('.o', '.f95', regex=False).tolist()
    return list_files


def merge_source(path_source_dir, path_target):
    """Short summary.

    Parameters
    ----------
    path_source_dir : string
        path to directory of dependencies
    path_target : string
        path for writing out the merged target file

    Returns
    -------
    path_target : Path
        path for writing out the merged target file

    """
    path_Makefile = os.path.join(path_source_dir, 'include.common')
    # get list of dependencies
    list_files = get_file_list(path_Makefile)

    f = open(path_target, 'w')
    for file in list_files:
        fp = open(os.path.join(path_source_dir, file), 'r')
        line = fp.readline()
        while line:
            # check if define wrf
            if line.lstrip().startswith('#ifdef wrf'):
                line = fp.readline()
                break_flag = False
                while break_flag == False:
                    if line.lstrip().startswith('#else'):
                        line = fp.readline()
                        while break_flag == False:
                            if line.lstrip().startswith('#endif'):
                                break_flag = True
                            else:
                                line = fp.readline()
                    elif line.lstrip().startswith('#endif'):
                        break_flag = True
                    else:
                        f.writelines(line)
                    line = fp.readline()
            # check if define nc
            elif line.lstrip().startswith('#ifdef nc'):
                line = fp.readline()
                break_flag = False
                while break_flag == False:
                    if line.lstrip().startswith('#else'):
                        line = fp.readline()
                        while break_flag == False:
                            if line.lstrip().startswith('#endif'):
                                break_flag = True
                            else:
                                f.writelines(line)
                                line = fp.readline()
                    elif line.lstrip().startswith('#endif'):
                        break_flag = True
                    line = fp.readline()
            else:
                f.writelines(line)
                line = fp.readline()
        fp.close()
        f.writelines('\n')
    f.close()

    return path_target


# path settings:
path_source_dir = '../SUEWS/SUEWS-SourceCode'
path_target = './module_sf_suewsdrv.F'


# merge files
# path_merged = merge_source(path_source_dir, path_target)
