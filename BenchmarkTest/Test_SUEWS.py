#!/usr/bin/env python

# from Benchmark_SUEWS import *
from shutil import copytree, rmtree, copyfile
import os
import sys
import numpy as np
import pandas as pd
from glob import glob
import f90nml
import errno
import filecmp
import tempfile
import itertools

# load global path
path_base=os.path.abspath(f90nml.read('BTS_config.nml')['file']
dir_input = path_base['dir_input'])
dir_exe = path_base['dir_exe'])
dir_baserun = path_base['dir_baserun'])

# load name of programme for testing
name_exe = f90nml.read('BTS_config.nml')['file']['name_exe']

# load physics options to test
dict_phy_opt_sel = f90nml.read('BTS_config.nml')['physics_test']

# %%auxiliary SUEWS functions
# suppress error info if needed:


class DevNull:
    def write(self, msg):
        pass


sys.stderr = DevNull()


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


# load results
def load_SUEWS_results(n_grid, n_year):
    path_out = os.path.join('Output/*SUEWS_60.txt')
    # re-order results into [year, grid] layout
    fl_res = np.array(
        sorted(glob(path_out))).reshape(
        n_grid, n_year).swapaxes(0, 1)
    res_sim0 = [np.array(
        [pd.read_csv(f_grid, sep='\s+', header=0).values
         for f_grid in fl_year]) for fl_year in fl_res]
    # re-order the results into [grid,time]
    res_sim = np.concatenate(res_sim0, axis=1)

    return res_sim


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


# geenerate a multi-grid-multi-year SiteSelect DataFrame
def gen_SiteSelect_multi(df_siteselect, n_grid):
    # process SiteSelect
    # we only test n_year=2 to see if multi-year run can be successful
    cfg_siteselect_base = df_siteselect.iloc[[0]]
    year0 = cfg_siteselect_base.Year[0]
    cfg_siteselect_base = pd.concat([cfg_siteselect_base, cfg_siteselect_base])
    cfg_siteselect_base.loc[:, 'Year'] = [year0, year0 + 1]
    cfg_siteselect_dim = cfg_siteselect_base.shape
    n_year = cfg_siteselect_dim[0]

    # create random surface fraction info for multiple grids
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
    df_siteselect_multi = pd.DataFrame(
        cfg_siteselect_x, columns=cfg_siteselect_base.columns)

    return df_siteselect_multi


# run simulation
def run_sim(name_sim, dir_exe, name_exe, dict_runcontrol, dict_initcond, df_siteselect,
            dir_save=tempfile.mkdtemp()):
    # TODO: support for user-specified forcing condition
    # create dir_save if not exisitng
    dir_save = os.path.abspath(os.path.expanduser(dir_save))
    if not os.path.exists(dir_save):
        os.mkdir(dir_save)
    # create folder `name_sim`
    dir_sys = os.getcwd()
    try:
        os.chdir(dir_save)
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
    yr_sim = df_siteselect.loc[:, 'Year'].unique().astype(int)
    grid_sim = df_siteselect.loc[:, 'Grid'].unique().astype(int)
    n_grid = len(grid_sim)
    n_year = len(yr_sim)

    # RunControl
    # set default values
    RunControl = {'runcontrol': dict_runcontrol.copy()}
    RunControl['runcontrol']['filecode'] = 'test'
    RunControl['runcontrol']['fileoutputpath'] = './Output/'

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
    # name_exe = 'SUEWS_V2018a'
    path_exe = os.path.join(dir_exe, name_exe)
    copyfile(path_exe, name_exe)
    os.chmod(name_exe, 755)

    # perform multi-grid run:
    # exit_code = os.system('./SUEWS_V2018a')
    # suppress output info
    os.system('./' + name_exe + ' &>/dev/null')

    # check if results generated:
    fl_output = glob('Output/*SUEWS_60.txt')
    if len(fl_output) > 0:
        # load results
        res_sim = load_SUEWS_results(n_grid, n_year)
        return res_sim
    else:
        # change back to original path
        os.chdir(dir_sys)
        return 'run failed!'

# %%

##############################################################################
# test functions for different purposes
# benchmark on the following tests:
# if compilation succeed?

# if running OK? i.e., generating expected files


# %%if single-grid multi-year run OK?
def test_multiyear(
        name_sim, dict_runcontrol, dict_initcond, df_siteselect,
        dir_save=tempfile.mkdtemp()):

    # generate a multi-grid SiteSelect file using only one grid
    df_siteselect_multi = gen_SiteSelect_multi(df_siteselect, n_grid=1)

    # process RunControl
    dict_runcontrol['resolutionfilesout'] = 3600
    dict_runcontrol['KeepTstepFilesOut'] = 0

    # result in [year, grid] order
    res_sim_multiyear = run_sim(
        name_sim, dir_exe, name_exe,
        dict_runcontrol, dict_initcond, df_siteselect_multi,
        dir_save)

    # test if multiple years have been run
    res_test = len(res_sim_multiyear.shape) > 0

    print('test_multiyear for',name_exe)
    print 'running here:', dir_save

    return res_test


# %%if multi-grid run produce the same resutls as their single-grid runs?
def test_multigrid(
        name_sim, dict_runcontrol, dict_initcond, df_siteselect,
        n_grid=3, dir_save=tempfile.mkdtemp()):
    # create folder `name_sim`
    dir_sys = os.getcwd()
    try:
        os.chdir(dir_save)
        # os.mkdir(name_sim)

        # os.mkdir('Output')
    except OSError as e:
        if e.errno == errno.EEXIST:
            print('Directory not created.')
        else:
            raise

    # generate a multi-grid SiteSelect file
    df_siteselect_multi = gen_SiteSelect_multi(df_siteselect, n_grid)

    # process RunControl
    dict_runcontrol['resolutionfilesout'] = 3600
    dict_runcontrol['KeepTstepFilesOut'] = 0

    # perform one multi-grid simulation and load results
    name_sim = 'multi_grid'
    # result in [year, grid] order
    res_sim_multigrid = run_sim(
        name_sim, dir_exe, name_exe,
        dict_runcontrol, dict_initcond,
        df_siteselect_multi,
        dir_save)

    # perform multiple single-grid simulation and load results
    res_sim_singlegrid = []
    list_grid = df_siteselect_multi.loc[:, 'Grid'].unique()
    for grid_name in list_grid:
        # grid_name = df_siteselect.loc[ind_grid, 'Grid']
        name_sim = 'single_grid_' + str(int(grid_name))
        res_sim_grid = run_sim(
            name_sim, dir_exe, name_exe,
            dict_runcontrol, dict_initcond,
            df_siteselect_multi.loc[df_siteselect_multi.Grid == grid_name],
            dir_save)
        res_sim_singlegrid.append(res_sim_grid)

    # combine `res_sim_singlegrid`
    res_sim_singlegrid = np.concatenate(tuple(res_sim_singlegrid))

    # test equality
    res_test = np.array_equal(res_sim_multigrid, res_sim_singlegrid)

    # change back to previous path
    os.chdir(dir_sys)

    print('test_multigrid for',name_exe)
    print 'running here:', dir_save


    return res_test


# %% if the test run can produce the same results as the SampleRun
def test_samerun(dir_baserun, dir_save=tempfile.mkdtemp()):
    dir_sys = os.getcwd()
    # load RunControl
    dict_runcontrol = load_SUEWS_RunControl(
        os.path.join(dir_baserun, 'RunControl.nml')
    )['runcontrol']
    # dir_input = dict_runcontrol['fileoutiutpath']
    dir_output = dict_runcontrol['fileoutputpath']

    # temp dir for testing
    dir_test = dir_save
    # clear tempdir
    rmtree(dir_test)

    # load configurations from dir_baserun
    copytree(dir_baserun, dir_test)
    os.chdir(dir_test)
    if os.path.exists(dir_output):
        rmtree(dir_output)
        os.mkdir(dir_output)

    # copy SUEWS executable
    os.remove(name_exe)
    path_exe = os.path.join(dir_exe, name_exe)
    copyfile(path_exe, name_exe)
    os.chmod(name_exe, 755)

    # perform simulation
    # exit_code = os.system('./SUEWS_V2018a')
    os.system('./' + name_exe + ' &>/dev/null')

    # compare results
    dir_res_sample = os.path.join(dir_baserun, dir_output)
    dir_res_test = os.path.join(dir_test, dir_output)
    # print('thse folders will be compared:')
    # print(dir_res_sample)
    # print(dir_res_test)

    common_files = [
        os.path.basename(x)
        for x in glob(os.path.join(dir_res_sample, '*'))]
    comp_files_test = filecmp.cmpfiles(
        dir_res_sample,
        dir_res_test,
        common_files,
        shallow=False)
    # print 'comp_files_test',comp_files_test
    # print 'common_files',common_files
    res_test = (set(comp_files_test[0]) == set(common_files))
    # print res_test
    if not res_test:
        # if not match, print mismatch list
        print 'these files are different:'
        print comp_files_test[1]


    dir_sys = os.chdir(dir_sys)
    # rmtree(dir_test)
    print('test_samerun for',name_exe)
    print 'running here:', dir_save

    return res_test


# %% diagnose working&faulty physics options
def test_physics(dict_runcontrol, dict_initcond, df_siteselect,
                 dict_phy_opt_sel,
                 dir_save=tempfile.mkdtemp()):

    print('test_physics for',name_exe)
    print 'running here:', dir_save

    # get options to test
    methods, options = zip(*dict_phy_opt_sel.items())
    options = [x if type(x) == list else [x] for x in options]
    list_to_test = [dict(zip(methods, v))
                    for v in itertools.product(*options)]

    # test selected physics schemes
    dict_test = {}
    for ind, cfg in enumerate(list_to_test):
        runcontrol_test = dict_runcontrol.copy()
        runcontrol_test.update(cfg)
        name_sim = str(ind)
        res_sim = run_sim(
            name_sim, dir_exe, name_exe,
            runcontrol_test,
            dict_initcond, df_siteselect,
            dir_save)
        dict_test.update({ind: res_sim})

    dict_test_OK = {k: 'fail' if type(v) == str else 'pass'
                    for k, v in dict_test.iteritems()}

    df_test = pd.DataFrame(list_to_test).assign(result=dict_test_OK.values())

    df_test.to_csv('~/Downloads/df_test.csv')
    # test results
    list_method_test = [c for c in df_test.columns if not c == 'result']
    df_test_pass = pd.concat(
        [df_test.loc[:, [c, 'result']].pivot_table(
            index='result', columns=c, aggfunc=len)
            for c in list_method_test],
        keys=list_method_test,
        axis=1).loc[
        'pass', :].to_frame().rename(
        columns={'pass': 'result'}).applymap(
        lambda x: 'pass' if x > 0 else 'fail')
    df_test_pass.index.set_names(['method', 'option'], inplace=True)

    # get `fail` options:
    list_fail = df_test_pass.loc[
        df_test_pass['result'] == 'fail'].index.tolist()

    return list_fail


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


##############################################################################
# def test_multigrid0(dir_sim, n_grid):
#     # run model with all grids
#     os.chdir(dir_sim)
#     if not os.path.exists('InputBase'):
#         copytree('Input', 'InputBase')
#
#     # process SiteSelect
#     fn_ss = 'SUEWS_SiteSelect.txt'
#     cfg_siteselect_base = pd.read_csv(
#         os.path.join('InputBase', fn_ss),
#         sep='\s+', skipfooter=2, header=1,
#         engine='python', comment='!')
#     cfg_siteselect_dim = cfg_siteselect_base.shape
#     n_year = cfg_siteselect_dim[0]
#
#     cfg_siteselect_header = np.vstack(
#         (np.arange(1, cfg_siteselect_dim[1] + 1), cfg_siteselect_base.columns))
#
#     # generate header for SiteSelect.txt
#     header_SS = '\n'.join(
#         [' '.join(line) for line in cfg_siteselect_header.astype(str)])
#
#     # create surface fraction info for multiple grids
#     # n_grid = 2  # number of grids
#     cfg_siteselect = np.tile(cfg_siteselect_base, (n_grid, 1))
#     fr_multigrid = [fr_line / np.sum(fr_line)
#                     for fr_line in np.random.rand(n_grid, 7)]
#
#     loc_fr = cfg_siteselect_base.columns.get_loc('Fr_Paved')
#     cfg_siteselect[:, loc_fr:loc_fr + 7] = np.repeat(
#         fr_multigrid, cfg_siteselect_dim[0], axis=0)
#
#     # specify grid numbers
#     cfg_siteselect[:, 0] = np.repeat(
#         np.arange(1, n_grid + 1), cfg_siteselect_dim[0])
#
#     # re-order results into [year, grid] layout
#     cfg_siteselect_x = cfg_siteselect.reshape(
#         n_grid, n_year, -1).swapaxes(0, 1).reshape(
#         n_grid * n_year, -1)
#
#     # create SiteSelect.txt
#     np.savetxt(os.path.join('Input', fn_ss), cfg_siteselect_x,
#                fmt=' '.join(['%i'] * 4 + ['%1.4f'] *
#                             (cfg_siteselect_dim[1] - 4)),
#                header=header_SS, footer='-9\n-9', comments='')
#
#     # clean output
#     if os.path.exists('Output'):
#         rmtree('Output')
#         os.mkdir('Output')
#
#     # perform multi-grid run:
#     exit_code = os.system('./SUEWS_V2018a')
#
#     # get results
#     # re-order results into [year, grid] layout
#     fl_res = np.array(
#         sorted(glob('Output/*SUEWS_60.txt'))).reshape(
#         n_grid, n_year).swapaxes(0, 1)
#     res_grid_multi0 = [np.array(
#         [pd.read_csv(f_grid, sep='\s+', header=0).values
#          for f_grid in fl_year]) for fl_year in fl_res]
#     # re-order the results into [grid,time]
#     res_grid_multi = np.concatenate(res_grid_multi0, axis=1)
#     res_grid_multi.shape
#     res_grid_multi[:, [0, -1], 0]
#     # perform single-grid run:
#     res_year_group = []
#     cfg_siteselect_year_grid = cfg_siteselect_x.reshape(n_year, n_grid, -1)
#     # cfg_siteselect_year_grid[0,:,0:2]
#     for cfg_year in cfg_siteselect_year_grid:
#         # save SiteSelect for individual grid
#         # cfg_year = cfg_siteselect_year_grid[0]
#         np.savetxt(os.path.join('Input', fn_ss), cfg_year,
#                    fmt=' '.join(['%i'] * 4 + ['%1.4f'] *
#                                 (cfg_siteselect_dim[1] - 4)),
#                    header=header_SS, footer='-9\n-9', comments='')
#         # clean output
#         if os.path.exists('Output'):
#             rmtree('Output')
#             os.mkdir('Output')
#         else:
#             os.mkdir('Output')
#         # return '0' means success
#         exit_code = os.system('./SUEWS_V2018a')
#         fl_res = sorted(glob('Output/*SUEWS_60.txt'))
#         res_year = np.array([pd.read_csv(
#             x, sep='\s+', header=0).values for x in fl_res])
#         res_year_group.append(res_year)
#
#     res_grid_single = np.concatenate(res_year_group, axis=1)
#     res_grid_single.shape
#     # recover the folder structure
#     if os.path.exists('Input'):
#         rmtree('Input')
#         copytree('InputBase', 'Input')
#     if os.path.exists('InputBase'):
#         rmtree('InputBase')
#         rmtree('Output')
#         os.mkdir('Output')
#
#     # test the results are equal
#     res_test = np.array_equal(res_grid_single, res_grid_multi)
#
#     return res_test






#
