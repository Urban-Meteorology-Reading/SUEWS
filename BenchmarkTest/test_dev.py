from Test_SUEWS import *
import unittest


# load basic configurations
# initial condition
dict_initcond = (load_SUEWS_nml(
    os.path.join(dir_input, 'InitialConditionstest_2004.nml')).to_dict()[
    'initialconditions'])
# runcontrol settings
dict_runcontrol = load_SUEWS_nml(
    os.path.join(dir_input, 'RunControl.nml')).to_dict()['runcontrol']
# siteselect info
df_siteselect = load_SUEWS_table(
    os.path.join(dir_input, 'SUEWS_SiteSelect.txt'))


# test case class for unit test
class Test_SUEWS(unittest.TestCase):
    def test_ok_multiyear(self):
        print('testing single-grid multi-year run ... ')
        name_sim = 'test-multi-year' + str(np.random.randint(10000))
        res_test = test_multiyear(
            name_sim, dict_runcontrol, dict_initcond, df_siteselect)
        self.assertTrue(res_test)

    def test_ok_multigrid(self):
        print('testing multi-grid multi-year run ... ')
        n_grid = 3
        name_sim = 'test-multi-grid' + str(np.random.randint(10000))
        res_test = test_multigrid(
            name_sim, dict_runcontrol, dict_initcond, df_siteselect,
            n_grid)
        self.assertTrue(res_test)

    def test_ok_samerun(self):
        print('testing if results could match the base run ... ')
        res_test = test_samerun(dir_baserun)
        self.assertTrue(res_test)


if __name__ == '__main__':
    unittest.main()
