# calibrate property values using optimisation
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
# from pandas import DataFrame as df
# import matplotlib.pyplot as plt


# construct an initilisation class
class spotpy_setup(object):
    def __init__(self):
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
        df_forcing_part = df_forcing.iloc[:288 * 30]
        dict_forcing_part = df_forcing_part.to_dict('index')

        # load observations
        list_file_Obs = glob.glob(os.path.join(
            dir_start, 'Obs/{}*txt'.format(filecode)))
        res_obs = pd.read_table(list_file_Obs[0])
        res_obs.rename(columns={'dectime': 'Dectime', 'QgHFP': 'QS_obs'},
                       inplace=True)

        self.res_obs = res_obs
        self.df_forcing_x = df_forcing_part
        self.dict_forcing_part = dict_forcing_part
        self.dict_state_init = dict_state_init
        # parameters to calibrate
        self.params = [
            spotpy.parameter.Uniform(
                'cpanohm', 0.1, 4),
            spotpy.parameter.Uniform(
                'chanohm', 0.1, 4),
            spotpy.parameter.Uniform(
                'kkanohm', 0.1, 2),
            spotpy.parameter.Uniform(
                'xwf', 0.1, 100)
        ]
        # self.params = [
        #     spotpy.parameter.Normal(
        #         'cpanohm', 1, 0.1, 0.01, 4),
        #     spotpy.parameter.Normal(
        #         'chanohm', 1.5, 0.35, 0.1, 4),
        #     spotpy.parameter.Normal(
        #         'kkanohm', 1, 0.3, 0.1, 2.0)
        # ]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # print 'vector',vector
        dict_state_init_x = self.dict_state_init[1].copy()
        dict_state_init_x['cpanohm'] = 1e6 * \
            vector[0] * np.ones((7), order='F')
        dict_state_init_x['chanohm'] = vector[1] * np.ones((7), order='F')
        dict_state_init_x['kkanohm'] = vector[2] * np.ones((7), order='F')
        dict_state_init_x['xwf'] = vector[3] * 1e-7

        dict_output, dict_state = sp.run_suews(
            self.dict_forcing_part, {1: dict_state_init_x})

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
        res_list_obs = res_comp.loc[:, 'QS_obs'].values.tolist()
        return res_list_obs  # should return a list

    def objectivefunction(self, simulation, evaluation):
        # print len(simulation)
        # print len(evaluation)
        # res_comp = evaluation.copy()
        # res_comp = res_comp.merge(simulation, on='Dectime')
        # res_list_obs = res_comp.loc[:, 'QS_obs'].values.tolist()
        # res_list_sim = res_comp.loc[:, 'QS_sim'].values.tolist()
        # print 'obs',res_list_obs
        # print 'sim',res_list_sim
        obj = - spotpy.objectivefunctions.mae(
            evaluation,
            simulation)
        return obj


reload(spotpy)
# sampling & calibration
xx_setup = spotpy_setup()
results = []
methods = [
    # 'mc',
    # 'mle',
    # 'lhs',
    'fscabc',  # good
    # 'sa',
    'sceua'
]
reps = 2
for method in methods:
    sampler = getattr(spotpy.algorithms,
                      method)(xx_setup,
                              dbname=method + '_SuPy',
                              dbformat='csv')
    try:
        sampler.sample(reps)
        results.append(sampler.getdata())
    except:
        pass

    # results.append(sampler.getdata())
#
# method = 'sceua'
# reps = 2
# sampler = getattr(spotpy.algorithms,
#                   method)(xx_setup,
#                           dbname=method + '_SuPy',
#                           dbformat='csv')
# # sampler.sample(reps, nChains=3)
# sampler.sample(reps)
# results.append(sampler.getdata())
# analysis
# spotpy.analyser.get_modelruns(results=results[0]);
# reload(spotpy)
spotpy.analyser.plot_bestmodelrun(results[0], xx_setup.evaluation())
# spotpy.analyser.plot_allmodelruns(results, xx_setup.evaluation())


spotpy.analyser.plot_regression(results[0], xx_setup.evaluation())


###############################################################################
# Calibrate water fluxes
#
###############################################################################

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
df_forcing_part = df_forcing.iloc[:288 * 30]
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
        self.df_forcing_x = df_forcing_x
        self.dict_forcing_part = df_forcing_x.to_dict('index')
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
                'xwf', 0.1, 100)
        ]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # print 'vector',vector
        dict_state_init_x = self.dict_state_init.copy()
        dict_state_init_x['cpanohm'] = 1e6 * \
            vector[0] * np.ones((7), order='F')
        dict_state_init_x['chanohm'] = vector[1] * np.ones((7), order='F')
        dict_state_init_x['kkanohm'] = vector[2] * np.ones((7), order='F')
        dict_state_init_x['xwf'] = vector[3] * 1e-7

        dict_output, dict_state = sp.run_suews(
            self.dict_forcing_part, {1: dict_state_init_x})

        # return dict_state
        self.dict_state = dict_state[1]

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
        res_list_obs = res_comp.loc[:, 'QS_obs'].values.tolist()
        return res_list_obs  # should return a list

    def objectivefunction(self, simulation, evaluation):
        # print len(simulation)
        # print len(evaluation)
        # res_comp = evaluation.copy()
        # res_comp = res_comp.merge(simulation, on='Dectime')
        # res_list_obs = res_comp.loc[:, 'QS_obs'].values.tolist()
        # res_list_sim = res_comp.loc[:, 'QS_sim'].values.tolist()
        # print 'obs',res_list_obs
        # print 'sim',res_list_sim
        obj = - spotpy.objectivefunctions.mae(
            evaluation,
            simulation)
        return obj


# day loop for calibration:
xdays = df_forcing_part.loc[:, 'id'].unique()
df_forcing_part_grp_list = df_forcing_part.groupby('id')
method = 'fscabc'
reps = 2
dict_state_init_day = dict_state_init[1]
for xday in xdays[:2]:
    df_forcing_day = df_forcing_part_grp_list.get_group(xday)

    xx_setup = spotpy_setup(dict_state_init=dict_state_init_day,
                            df_forcing_x=df_forcing_day,
                            df_obs=res_obs)
    sampler = getattr(spotpy.algorithms,
                      method)(xx_setup,
                              dbname=method + '_SuPy',
                              dbformat='csv')
    try:
        sampler.sample(reps)
        # results.append(sampler.getdata())
        dict_state_init_day = sampler.dict_state.values()[-1]
