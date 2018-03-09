###############################################################################
# Calibrate water fluxes
#
###############################################################################
# Load the SPOT package into your working storage
import spotpy
# Load the Plotting extension
from spotpy import analyser
# Import the two dimensional Rosenbrock example
# from spotpy.examples.spot_setup_rosenbrock import spot_setup
# Import SuPy module: SUEWS model
# load dependency modules
import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
import copy
# from pandas import DataFrame as df
# import matplotlib.pyplot as plt

# initialise SUEWS settings
# dir_input = './Input'
dir_start = '../SampleRun'
dict_mod_cfg, dict_state_init = sp.init_SUEWS_dict(dir_start)

# load met forcing
filecode = dict_mod_cfg['filecode']
tstep = dict_state_init[1]['tstep']
list_file_MetForcing = glob.glob(os.path.join(
    dir_start, 'Input/{}*{}*txt'.format(filecode, tstep / 60)))
# load as DataFrame:
df_forcing = sp.load_SUEWS_MetForcing_df(list_file_MetForcing[1])

# prepare forcing:
df_forcing_part = df_forcing.iloc[:]
dict_forcing_part = df_forcing_part.to_dict('index')

# load observations
list_file_Obs = glob.glob(os.path.join(
    dir_start, 'Obs/{}*txt'.format(filecode)))
res_obs = pd.read_table(list_file_Obs[0])
res_obs.rename(columns={'dectime': 'Dectime', 'QgHFP': 'QS_obs'}, inplace=True)


# construct an initilisation class
class spotpy_setup(object):
    def __init__(self, dict_state_init, df_forcing_x, df_obs):
        # initialise SUEWS settings
        self.dict_state_init = dict_state_init

        # forcing passed as initial setting
        self.df_forcing_x = df_forcing_x.copy()
        # print 'self.df_forcing_x.shape',self.df_forcing_x.shape
        # self.dict_forcing_part = self.df_forcing_x.to_dict('index')
        # print 'dict_forcing_part.keys',self.dict_forcing_part.keys()
        # observations passed as initial settings
        self.res_obs = df_obs

        # parameters to calibrate
        self.params = [
            spotpy.parameter.Uniform(
                'cpanohm', 0.1, 4),
            spotpy.parameter.Uniform(
                'chanohm', 0.1, 4),
            spotpy.parameter.Uniform(
                'kkanohm', 0.1, 2),
            spotpy.parameter.Uniform(
                'xwf', -100, 100)
        ]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # print 'vector', vector
        dict_state_init_x = self.dict_state_init.copy()
        dict_state_init_x['cpanohm'] = 1e6 * \
            vector[0] * np.ones((7), order='F')
        dict_state_init_x['chanohm'] = vector[1] * np.ones((7), order='F')
        dict_state_init_x['kkanohm'] = vector[2] * np.ones((7), order='F')
        dict_state_init_x['xwf'] = vector[3] * 1e-7

        # xx = np.array([dict_state_init_x[k]
        #                for k in ['cpanohm', 'kkanohm', 'chanohm', 'xwf']])
        # print(np.array(xx,dtype=np.float).shape)
        # print sorted(dict_state_init_x.keys())
        # print 'properties before pass-in'
        # print dict_state_init_x['cpanohm'][0]
        # print dict_state_init_x['chanohm'][0]
        # print dict_state_init_x['kkanohm'][0]
        # print dict_state_init_x['xwf']

        # force change `id` at each simulation
        # to trigger AnOHM recalculate coefficients
        # df_forcing_x = self.df_forcing_x

        # print 'before'
        # print 'day0', self.df_forcing_x.iloc[0]['id']
        df_forcing_x = self.df_forcing_x
        # xday = df_forcing_x.iloc[0]['id'] + 3
        # # print 'xday', xday
        # df_forcing_x['id'] = xday
        # met_grid = df_forcing_x.iloc[0]['metforcingdata_grid']
        # met_grid['id'] = xday
        # df_forcing_x.loc[:, 'metforcingdata_grid'] = met_grid
        # df_forcing_x.loc[:, 'metforcingdata_grid'] = df_forcing_x.loc[
        #     :, 'metforcingdata_grid'].map(lambda x: met_grid)

        # print 'after'
        # print df_forcing_x.iloc[1]['id']
        # print df_forcing_x.iloc[1]['metforcingdata_grid'].iloc[1]['id']

        dict_output, dict_state = sp.run_suews(
            df_forcing_x, {1: dict_state_init_x})

        # return dict_state
        self.dict_state = dict_state

        df_output = sp.pack_df_output(dict_output)
        res_sim = df_output.loc[1, 'SUEWS'].loc[:, ['Dectime', 'QS']]
        res_sim.rename(columns={'QS': 'QS_sim'}, inplace=True)

        res_comp = self.res_obs
        res_comp = res_comp.merge(res_sim, on='Dectime')
        # res_list_obs = res_comp.loc[:, 'QS_obs'].values.tolist()
        res_list_sim = res_comp.loc[:, 'QS_sim'].values.tolist()

        return res_list_sim  # should return a list

    def evaluation(self):
        res_comp = self.df_forcing_x.copy()
        res_comp.rename(columns={'dectime': 'Dectime'}, inplace=True)
        res_comp = res_comp.merge(self.res_obs, on='Dectime')
        # print res_comp.loc[:,['Dectime','QS_sim','QS_obs']]
        res_list_obs = res_comp.loc[:, 'QS_obs'].values.tolist()
        return res_list_obs  # should return a list

    def objectivefunction(self, simulation, evaluation):
        # print len(simulation)
        # print len(evaluation)
        # res_comp = evaluation.copy()
        # res_comp = res_comp.merge(simulation, on='Dectime')
        # res_list_obs = res_comp.loc[:, 'QS_obs'].values.tolist()
        # res_list_sim = res_comp.loc[:, 'QS_sim'].values.tolist()
        # print 'obs', evaluation
        # print 'sim', simulation
        obj = - spotpy.objectivefunctions.mae(
            evaluation,
            simulation)
        # obj = spotpy.objectivefunctions.nashsutcliffe(
        #     evaluation,
        #     simulation)
        # print 'obj', obj
        # print '\n'
        return obj


# day loop for calibration:
xdays = df_forcing_part.loc[:, 'id'].unique()
df_forcing_part_grp_list = df_forcing_part.groupby('id')
method = 'sceua'
reps = 40
dict_state_init_day = dict_state_init[1]
for xday in xdays[:90]:
    print 'xday', xday
    # print 'dict_state_init_day',dict_state_init_day.keys()
    df_forcing_day = df_forcing_part_grp_list.get_group(xday)

    xx_setup = spotpy_setup(dict_state_init=copy.deepcopy(dict_state_init_day),
                            df_forcing_x=df_forcing_day.copy(),
                            df_obs=res_obs)
    sampler = getattr(spotpy.algorithms,
                      method)(xx_setup,
                              dbname='{method}_SuPy_{id}'.format(
                                  method=method, id=xday),
                              dbformat='csv')
    # try:
    sampler.sample(reps)
    t_last = sorted(xx_setup.dict_state.keys())[-1]
    print 't_last', t_last
    dict_state_init_day = copy.deepcopy(xx_setup.dict_state[t_last][1])
    spotpy.analyser.plot_bestmodelrun(sampler.getdata(), xx_setup.evaluation())
    spotpy.analyser.plot_bestmodelrun(
        sampler.getdata(), xx_setup.evaluation())
    # except:
    #     pass

# sampler.get_parameters()
# spotpy.analyser.plot_bestmodelrun(sampler.getdata(), xx_setup.evaluation())
