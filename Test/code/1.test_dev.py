from shutil import rmtree
import Test_SUEWS as ts
import f90nml
import os
import numpy as np
import unittest
from pathlib import Path
from tempfile import gettempdir, TemporaryDirectory
from shutil import copyfile, copytree

fn_nml = 'BTS_config.nml'
# load basic configurations
# load path
nml = f90nml.read(fn_nml)
cfg_file = nml['file']
dir_exe, path_baserun = (
    os.path.abspath(cfg_file[x])
    for x in ['dir_exe', 'dir_baserun'])
path_baserun = Path(path_baserun)
# dir_exe = os.path.abspath(path_base[])
# dir_baserun = path_base[]


# copy all input files to the release folder
path_release_input = Path('../../Release/InputTables').resolve()
# version specific input folder under release folder
path_input_ver = (path_release_input / path_baserun.name).resolve()
print('dir_input:', path_input_ver)
# clean existing files
if path_input_ver.exists():
    print('cleaning existing files...')
    rmtree(path_input_ver)
# make an empty directory
path_input_ver.mkdir()


# copy runcontrol
path_runctrl_base = path_baserun/'RunControl.nml'
print('path_runctrl_base:', path_runctrl_base)
path_runctrl_input = path_input_ver/'RunControl.nml'
print('path_runctrl_input:', path_runctrl_input)
copyfile(path_runctrl_base, path_runctrl_input)
dict_runcontrol = ts.load_SUEWS_nml(path_runctrl_base)['runcontrol']
# copy other input tables and initial conditions
path_base_input = (path_baserun / dict_runcontrol['fileinputpath'])
for x in path_base_input.glob('*'):
    if x.is_dir():
        copytree(x, path_input_ver / x.name)
    else:
        copyfile(x, path_input_ver / x.name)

# load name of programme for testing
name_exe = cfg_file['name_exe']

# load test configurations
cfg_test = nml['test']
flag_multi_grid = True if cfg_test['multi_grid'] == 1 else False
flag_multi_year = True if cfg_test['multi_year'] == 1 else False
flag_same_run = True if cfg_test['same_run'] == 1 else False
flag_test_phys = True if cfg_test['test_phys'] == 1 else False
flag_test_complete = True if cfg_test['test_complete'] == 1 else False
test_number = cfg_test['test_number']

# load physics options to test
dict_phy_opt_sel = nml['physics_test']

# runcontrol settings
dict_runcontrol = ts.load_SUEWS_nml(
    os.path.join(path_input_ver, 'RunControl.nml')).to_dict()['runcontrol']
# initial condition
dict_initcond = (ts.load_SUEWS_nml(
    os.path.join(path_input_ver, 'InitialConditionstest_2004.nml')).to_dict()[
    'initialconditions'])
# siteselect info
df_siteselect = ts.load_SUEWS_table(
    os.path.join(path_input_ver, 'SUEWS_SiteSelect.txt'))


# test case class for unit test
class Test_SUEWS(unittest.TestCase):
    def test_ok_multiyear(self):
        print('***************************************')
        if flag_multi_year:
            print('testing single-grid multi-year run ... ')
            name_sim = 'test-multi-year' + str(np.random.randint(10000))
            res_test = ts.test_multiyear(
                name_sim, name_exe, dict_runcontrol, dict_initcond, df_siteselect,
                dir_exe, path_input_ver)
            self.assertTrue(res_test)
        else:
            print('single-grid multi-year test skipped ... ')
        print('  ')
        print('***************************************')

    def test_ok_multigrid(self):
        print('')
        print('***************************************')
        if flag_multi_grid:
            print('testing multi-grid multi-year run ... ')
            n_grid = 3
            name_sim = 'test-multi-grid' + str(np.random.randint(10000))
            res_test = ts.test_multigrid(
                name_sim, name_exe,
                dict_runcontrol, dict_initcond, df_siteselect,
                n_grid, dir_exe, path_input_ver)
            self.assertTrue(res_test)
        else:
            print('multi-grid multi-year test skipped ... ')
        print('  ')
        print('***************************************')

    def test_ok_samerun(self):
        print('')
        print('****************************************************')

        if flag_same_run:
            print('testing if results could match the standard run ... ')
            print('N.B.: THE TEST IS GENERATED USING LONGTERM FORCING!')
            print('N.B.: test_2004_data_60.txt.long')
            name_sim = 'test-same-run' + str(np.random.randint(10000))
            res_test = ts.test_samerun(name_sim, name_exe,
                                        dict_runcontrol, dict_initcond,
                                        df_siteselect,
                                        dir_exe, path_baserun)
            self.assertTrue(res_test)
        else:
            print('identity test skipped ... ')
        print('  ')
        print('****************************************************')

    def test_ok_physics(self):
        print('')
        print('************************************************')
        if flag_test_phys:
            print('testing if some physics schemes are working ... ')
            if flag_test_complete:
                if test_number>0:
                    print(
                        f'testing in selective mode: a randomly chosen {test_number} combinations of physics schemes will be tested!')
                else:
                    print('testing in complete mode: all physics schemes will be tested!')
            else:
                print('testing in concise mode: only part of physics schemes will be tested!')
            # show options to test
            res_list_fail = ts.test_physics(
                name_exe, path_input_ver, dir_exe,
                dict_runcontrol, dict_initcond, df_siteselect,
                dict_phy_opt_sel,
                flag_test_complete,
                test_number
            )

            # `0` means no failure: all options can pass test
            res_test = len(res_list_fail) == 0
            print(res_list_fail)
            # print out all faulty options
            if not res_test:
                print('faulty options found:')
                for fail in res_list_fail:
                    print(fail)

            self.assertTrue(res_test)
        else:
            print('physics test skipped ... ')
        print('  ')
        print('************************************************')


if __name__ == '__main__':
    unittest.main()
