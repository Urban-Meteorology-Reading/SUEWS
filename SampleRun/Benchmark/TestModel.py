#!/usr/bin/env python

# from Benchmark_SUEWS import *
from shutil import copytree, rmtree, copyfile
import os
import numpy as np
import pandas as pd
from glob import glob
import f90nml
import errno


dir_base = '/Users/sunt05/Dropbox/8-Research/98.ReadingWork/10.SUEWS-FORTRAN/'
dir_input = os.path.join(dir_base, 'SampleRun/Benchmark/BaseInput')
dir_exe = os.path.join(dir_base, 'ReleaseRepo/build/macOS')


# benchmark on the following tests:
# if compilation succeed?

# if running OK? i.e., generating expected files

# load model settings
# load configurations: mod_config
# process RunControl.nml
# this function can handle all SUEWS nml files


def load_SUEWS_nml(xfile):
    # remove case issues
    xfile = path_insensitive(xfile)
    df = pd.DataFrame(f90nml.read(xfile))
    return df


def load_SUEWS_RunControl(xfile):
    lib_RunControl = load_SUEWS_nml(xfile)
    # return DataFrame containing settings
    return lib_RunControl


# load all tables (xgrid.e., txt files)
def load_SUEWS_table(fileX):
    # remove case issues
    fileX = path_insensitive(fileX)
    rawdata = pd.read_table(fileX, delim_whitespace=True,
                            comment='!', error_bad_lines=True,
                            skiprows=1).dropna()
    return rawdata


# save SiteSelect to file
def save_SiteSelect(df_siteselect, fn_ss='Input/SUEWS_SiteSelect.txt'):
    # process SiteSelect
    # fn_ss = 'SUEWS_SiteSelect.txt'
    cfg_siteselect_base = df_siteselect
    cfg_siteselect_dim = cfg_siteselect_base.shape
    cfg_siteselect_x = df_siteselect.values
    # n_year = cfg_siteselect_dim[0]

    cfg_siteselect_header = np.vstack(
        (np.arange(1, cfg_siteselect_dim[1] + 1), cfg_siteselect_base.columns))

    # generate header for SiteSelect.txt
    header_SS = '\n'.join(
        [' '.join(line) for line in cfg_siteselect_header.astype(str)])

    # create SiteSelect.txt
    np.savetxt(fn_ss, cfg_siteselect_x,
               fmt=' '.join(['%i'] * 4 + ['%1.4f'] *
                            (cfg_siteselect_dim[1] - 4)),
               header=header_SS, footer='-9\n-9', comments='')
    # print('SiteSelect saved!')


# save InitCond to file
def save_InitCond(dict_initcond, year, grid=''):
    InitCond = {'initialconditions': dict_initcond.copy()}
    for k, v in InitCond['initialconditions'].iteritems():
        if int(v) == v:
            InitCond['initialconditions'][k] = int(v)

    # save RunControl to file
    fn_nml = os.path.join(
        'Input',
        'InitialConditionstest{grid}_{year}.nml').format(
        grid=grid,
        year=year)
    if os.path.exists(fn_nml):
        os.remove(fn_nml)
    f90nml.write(InitCond, fn_nml)

    print('InitCond saved!')


# run simulation
def run_sim(name_sim, dict_runcontrol, dict_initcond, df_siteselect):
    # TODO: support for user-specified forcing condition
    # create folder `name_sim`
    dir_save = os.getcwd()
    try:
        os.mkdir(name_sim)
        os.chdir(name_sim)
        os.mkdir('Output')
    except OSError as e:
        if e.errno == errno.EEXIST:
            print('Directory not created.')
        else:
            raise
    # copy base input files
    copytree(dir_input, 'Input')

    # get sim info
    yr_sim = df_siteselect.Year.unique().astype(int)
    grid_sim = df_siteselect.Grid.unique().astype(int)
    n_grid = len(grid_sim)
    n_year = len(yr_sim)

    # RunControl
    # set default values
    RunControl = {'runcontrol': dict_runcontrol.copy()}
    RunControl['runcontrol']['filecode'] = 'test'
    RunControl['runcontrol']['fileoutputpath'] = './Output'

    # Initial condition
    # multi-grid initcond or not?
    if np.array_equal(grid_sim, dict_initcond.keys()):
        # use specific InitCond
        RunControl['runcontrol']['multipleinitfiles'] = 1
        # set default values
        for grid in grid_sim:
            save_InitCond(dict_initcond[grid], min(yr_sim), grid)
    else:
        # use same InitCond
        RunControl['runcontrol']['multipleinitfiles'] = 0
        save_InitCond(dict_initcond, min(yr_sim))

    # save SiteSelect to file
    save_SiteSelect(df_siteselect)

    # save RunControl to file
    f90nml.write(RunControl, 'RunControl.nml')

    # copy SUEWS executable
    name_exe = 'SUEWS_V2018a'
    copyfile(os.path.join(dir_exe, name_exe), name_exe)
    os.chmod(name_exe, 755)

    # perform multi-grid run:
    # exit_code = os.system('./SUEWS_V2018a')
    os.system('./SUEWS_V2018a')

    # load results
    # re-order results into [year, grid] layout
    fl_res = np.array(
        sorted(glob('Output/*SUEWS_60.txt'))).reshape(
        n_grid, n_year).swapaxes(0, 1)
    res_sim0 = [np.array(
        [pd.read_csv(f_grid, sep='\s+', header=0).values
         for f_grid in fl_year]) for fl_year in fl_res]
    # re-order the results into [grid,time]
    res_sim = np.concatenate(res_sim0, axis=1)

    # change back to original path
    os.chdir(dir_save)

    return res_sim


# if single-grid multi-year run OK?

# if multi-grid run produce the same resutls as their single-grid runs?
def test_multigrid(
        name_sim, dict_runcontrol, dict_initcond, df_siteselect,
        n_grid):
    # create folder `name_sim`
    dir_save = os.getcwd()
    try:
        os.mkdir(name_sim)
        os.chdir(name_sim)
        os.mkdir('Output')
    except OSError as e:
        if e.errno == errno.EEXIST:
            print('Directory not created.')
        else:
            raise

    # process SiteSelect
    cfg_siteselect_base = df_siteselect.iloc[[0]]
    year0 = cfg_siteselect_base.Year[0]
    cfg_siteselect_base = pd.concat([cfg_siteselect_base, cfg_siteselect_base])
    cfg_siteselect_base.loc[:, 'Year'] = [year0, year0 + 1]
    cfg_siteselect_dim = cfg_siteselect_base.shape
    n_year = cfg_siteselect_dim[0]

    # create surface fraction info for multiple grids
    # n_grid = 2  # number of grids
    cfg_siteselect = np.tile(cfg_siteselect_base.values, (n_grid, 1))
    fr_multigrid = [fr_line / np.sum(fr_line)
                    for fr_line in np.random.rand(n_grid, 7)]

    loc_fr = cfg_siteselect_base.columns.get_loc('Fr_Paved')
    cfg_siteselect[:, loc_fr:loc_fr + 7] = np.repeat(
        fr_multigrid, cfg_siteselect_dim[0], axis=0)

    # specify grid numbers
    cfg_siteselect[:, 0] = np.repeat(
        np.arange(1, n_grid + 1), cfg_siteselect_dim[0])

    # re-order results into [year, grid] layout
    cfg_siteselect_x = cfg_siteselect.reshape(
        n_grid, n_year, -1).swapaxes(0, 1).reshape(
        n_grid * n_year, -1)

    # perform multi-grid run:
    df_siteselect = pd.DataFrame(
        cfg_siteselect_x, columns=cfg_siteselect_base.columns)

    # process InitCond
    # dict_initcond = (load_SUEWS_nml(
    #     'Benchmark/BaseInput/InitialConditionstest_2004.nml').to_dict()[
    #     'initialconditions'])

    # process RunControl
    dict_runcontrol['resolutionfilesout'] = 3600
    dict_runcontrol['KeepTstepFilesOut'] = 0

    # perform simulation
    name_sim = 'multi-grid'
    exit_code = run_sim(
        name_sim, dict_runcontrol, dict_initcond, df_siteselect)

    # # get results
    # # re-order results into [year, grid] layout
    # fl_res = np.array(
    #     sorted(glob('Output/*SUEWS_60.txt'))).reshape(
    #     n_grid, n_year).swapaxes(0, 1)
    # res_grid_multi0 = [np.array(
    #     [pd.read_csv(f_grid, sep='\s+', header=0).values
    #      for f_grid in fl_year]) for fl_year in fl_res]
    # # re-order the results into [grid,time]
    # res_grid_multi = np.concatenate(res_grid_multi0, axis=1)
    # res_grid_multi.shape
    # res_grid_multi[:, [0, -1], 0]
    # # perform single-grid run:
    # res_year_group = []
    # cfg_siteselect_year_grid = cfg_siteselect_x.reshape(n_year, n_grid, -1)
    # # cfg_siteselect_year_grid[0,:,0:2]
    # for cfg_year in cfg_siteselect_year_grid:
    #     # save SiteSelect for individual grid
    #     # cfg_year = cfg_siteselect_year_grid[0]
    #     np.savetxt(os.path.join('Input', fn_ss), cfg_year,
    #                fmt=' '.join(['%i'] * 4 + ['%1.4f'] *
    #                             (cfg_siteselect_dim[1] - 4)),
    #                header=header_SS, footer='-9\n-9', comments='')
    #     # clean output
    #     if os.path.exists('Output'):
    #         rmtree('Output')
    #         os.mkdir('Output')
    #     else:
    #         os.mkdir('Output')
    #     # return '0' means success
    #     exit_code = os.system('./SUEWS_V2018a')
    #     fl_res = sorted(glob('Output/*SUEWS_60.txt'))
    #     res_year = np.array([pd.read_csv(
    #         x, sep='\s+', header=0).values for x in fl_res])
    #     res_year_group.append(res_year)
    #
    # res_grid_single = np.concatenate(res_year_group, axis=1)
    # res_grid_single.shape
    # # recover the folder structure
    # if os.path.exists('Input'):
    #     rmtree('Input')
    #     copytree('InputBase', 'Input')
    # if os.path.exists('InputBase'):
    #     rmtree('InputBase')
    #     rmtree('Output')
    #     os.mkdir('Output')
    #
    # # test the results are equal
    # res_test = np.array_equal(res_grid_single, res_grid_multi)

    return exit_code


def test_multigrid0(dir_sim, n_grid):
    # run model with all grids
    os.chdir(dir_sim)
    if not os.path.exists('InputBase'):
        copytree('Input', 'InputBase')

    # process SiteSelect
    fn_ss = 'SUEWS_SiteSelect.txt'
    cfg_siteselect_base = pd.read_csv(
        os.path.join('InputBase', fn_ss),
        sep='\s+', skipfooter=2, header=1,
        engine='python', comment='!')
    cfg_siteselect_dim = cfg_siteselect_base.shape
    n_year = cfg_siteselect_dim[0]

    cfg_siteselect_header = np.vstack(
        (np.arange(1, cfg_siteselect_dim[1] + 1), cfg_siteselect_base.columns))

    # generate header for SiteSelect.txt
    header_SS = '\n'.join(
        [' '.join(line) for line in cfg_siteselect_header.astype(str)])

    # create surface fraction info for multiple grids
    # n_grid = 2  # number of grids
    cfg_siteselect = np.tile(cfg_siteselect_base, (n_grid, 1))
    fr_multigrid = [fr_line / np.sum(fr_line)
                    for fr_line in np.random.rand(n_grid, 7)]

    loc_fr = cfg_siteselect_base.columns.get_loc('Fr_Paved')
    cfg_siteselect[:, loc_fr:loc_fr + 7] = np.repeat(
        fr_multigrid, cfg_siteselect_dim[0], axis=0)

    # specify grid numbers
    cfg_siteselect[:, 0] = np.repeat(
        np.arange(1, n_grid + 1), cfg_siteselect_dim[0])

    # re-order results into [year, grid] layout
    cfg_siteselect_x = cfg_siteselect.reshape(
        n_grid, n_year, -1).swapaxes(0, 1).reshape(
        n_grid * n_year, -1)

    # create SiteSelect.txt
    np.savetxt(os.path.join('Input', fn_ss), cfg_siteselect_x,
               fmt=' '.join(['%i'] * 4 + ['%1.4f'] *
                            (cfg_siteselect_dim[1] - 4)),
               header=header_SS, footer='-9\n-9', comments='')

    # clean output
    if os.path.exists('Output'):
        rmtree('Output')
        os.mkdir('Output')

    # perform multi-grid run:
    exit_code = os.system('./SUEWS_V2018a')

    # get results
    # re-order results into [year, grid] layout
    fl_res = np.array(
        sorted(glob('Output/*SUEWS_60.txt'))).reshape(
        n_grid, n_year).swapaxes(0, 1)
    res_grid_multi0 = [np.array(
        [pd.read_csv(f_grid, sep='\s+', header=0).values
         for f_grid in fl_year]) for fl_year in fl_res]
    # re-order the results into [grid,time]
    res_grid_multi = np.concatenate(res_grid_multi0, axis=1)
    res_grid_multi.shape
    res_grid_multi[:, [0, -1], 0]
    # perform single-grid run:
    res_year_group = []
    cfg_siteselect_year_grid = cfg_siteselect_x.reshape(n_year, n_grid, -1)
    # cfg_siteselect_year_grid[0,:,0:2]
    for cfg_year in cfg_siteselect_year_grid:
        # save SiteSelect for individual grid
        # cfg_year = cfg_siteselect_year_grid[0]
        np.savetxt(os.path.join('Input', fn_ss), cfg_year,
                   fmt=' '.join(['%i'] * 4 + ['%1.4f'] *
                                (cfg_siteselect_dim[1] - 4)),
                   header=header_SS, footer='-9\n-9', comments='')
        # clean output
        if os.path.exists('Output'):
            rmtree('Output')
            os.mkdir('Output')
        else:
            os.mkdir('Output')
        # return '0' means success
        exit_code = os.system('./SUEWS_V2018a')
        fl_res = sorted(glob('Output/*SUEWS_60.txt'))
        res_year = np.array([pd.read_csv(
            x, sep='\s+', header=0).values for x in fl_res])
        res_year_group.append(res_year)

    res_grid_single = np.concatenate(res_year_group, axis=1)
    res_grid_single.shape
    # recover the folder structure
    if os.path.exists('Input'):
        rmtree('Input')
        copytree('InputBase', 'Input')
    if os.path.exists('InputBase'):
        rmtree('InputBase')
        rmtree('Output')
        os.mkdir('Output')

    # test the results are equal
    res_test = np.array_equal(res_grid_single, res_grid_multi)

    return res_test


##############################################################################
# auxiliary functions
# resolve path case issues
def path_insensitive(path):
    """
    Get a case-insensitive path for use on a case sensitive system.

    >>> path_insensitive('/Home')
    '/home'
    >>> path_insensitive('/Home/chris')
    '/home/chris'
    >>> path_insensitive('/HoME/CHris/')
    '/home/chris/'
    >>> path_insensitive('/home/CHRIS')
    '/home/chris'
    >>> path_insensitive('/Home/CHRIS/.gtk-bookmarks')
    '/home/chris/.gtk-bookmarks'
    >>> path_insensitive('/home/chris/.GTK-bookmarks')
    '/home/chris/.gtk-bookmarks'
    >>> path_insensitive('/HOME/Chris/.GTK-bookmarks')
    '/home/chris/.gtk-bookmarks'
    >>> path_insensitive("/HOME/Chris/I HOPE this doesn't exist")
    "/HOME/Chris/I HOPE this doesn't exist"
    """

    return _path_insensitive(path) or path


def _path_insensitive(path):
    """
    Recursive part of path_insensitive to do the work.
    """

    if path == '' or os.path.exists(path):
        return path

    base = os.path.basename(path)  # may be a directory or a file
    dirname = os.path.dirname(path)

    suffix = ''
    if not base:  # dir ends with a slash?
        if len(dirname) < len(path):
            suffix = path[:len(path) - len(dirname)]

        base = os.path.basename(dirname)
        dirname = os.path.dirname(dirname)

    if not os.path.exists(dirname):
        dirname = _path_insensitive(dirname)
        if not dirname:
            return

    # at this point, the directory exists but not the file

    try:  # we are expecting dirname to be a directory, but it could be a file
        files = os.listdir(dirname)
    except OSError:
        return

    baselow = base.lower()
    try:
        basefinal = next(fl for fl in files if fl.lower() == baselow)
    except StopIteration:
        return

    if basefinal:
        return os.path.join(dirname, basefinal) + suffix
    else:
        return
