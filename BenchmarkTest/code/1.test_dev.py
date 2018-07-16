import Test_SUEWS as ts
import f90nml
import os
import numpy as np
import unittest

fn_nml = 'BTS_config.nml'
# load basic configurations
# load path
nml = f90nml.read(fn_nml)
cfg_file = nml['file']
dir_input, dir_exe, dir_baserun = (
    os.path.abspath(cfg_file[x])
    for x in ['dir_input', 'dir_exe', 'dir_baserun'])
# dir_exe = os.path.abspath(path_base[])
# dir_baserun = path_base[]

# load name of programme for testing
name_exe = cfg_file['name_exe']

# load physics options to test
dict_phy_opt_sel = nml['physics_test']

# runcontrol settings
dict_runcontrol = ts.load_SUEWS_nml(
    os.path.join(dir_input, 'RunControl.nml')).to_dict()['runcontrol']
# initial condition
dict_initcond = (ts.load_SUEWS_nml(
    os.path.join(dir_input, 'InitialConditionstest_2004.nml')).to_dict()[
    'initialconditions'])
# siteselect info
df_siteselect = ts.load_SUEWS_table(
    os.path.join(dir_input, 'SUEWS_SiteSelect.txt'))


# test case class for unit test
class Test_SUEWS(unittest.TestCase):
    def test_ok_multiyear(self):
        print('***************************************')
        print('testing single-grid multi-year run ... ')
        name_sim = 'test-multi-year' + str(np.random.randint(10000))
        res_test = ts.test_multiyear(
            name_sim, name_exe, dict_runcontrol, dict_initcond, df_siteselect,
            dir_exe, dir_input)
        self.assertTrue(res_test)
        print('  ')
        print('***************************************')

    def test_ok_multigrid(self):
        print('')
        print('***************************************')
        print('testing multi-grid multi-year run ... ')
        n_grid = 3
        name_sim = 'test-multi-grid' + str(np.random.randint(10000))
        res_test = ts.test_multigrid(
            name_sim, name_exe,
            dict_runcontrol, dict_initcond, df_siteselect,
            n_grid, dir_exe, dir_input)
        self.assertTrue(res_test)
        print('  ')
        print('***************************************')

    def test_ok_samerun(self):
        print('')
        print('****************************************************')
        print('testing if results could match the standard run ... ')
        name_sim = 'test-same-run' + str(np.random.randint(10000))
        res_test = ts.test_samerun(name_sim, name_exe,
                                   dict_runcontrol, dict_initcond,
                                   df_siteselect,
                                   dir_exe, dir_baserun)
        self.assertTrue(res_test)
        print('  ')
        print('****************************************************')

    def test_ok_physics(self):
        print('')
        print('************************************************')
        print('testing if some physics schemes are working ... ')
        # show options to test
        res_list_fail = ts.test_physics(
            name_exe, dir_input, dir_exe,
            dict_runcontrol, dict_initcond, df_siteselect,
            dict_phy_opt_sel)

        # `0` means no failure: all options can pass test
        res_test = len(res_list_fail) == 0
        print(res_list_fail)
        # print out all faulty options
        if not res_test:
            print('faulty options found:')
            for fail in res_list_fail:
                print(fail)

        self.assertTrue(res_test)
        print('  ')
        print('************************************************')


if __name__ == '__main__':
    unittest.main()
