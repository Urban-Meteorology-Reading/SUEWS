import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib.pyplot as plt
import copy
# import collections
reload(sp)

# initialise SUEWS settings
# dir_input = './Input'
dir_start = '../SampleRun'
dict_mod_cfg, dict_state_init = sp.init_SUEWS_dict(dir_start)

# load met forcing
filecode = dict_mod_cfg['filecode']
tstep = dict_state_init[1]['tstep']
metfile_pattern='Input/{}*{}*txt'.format(filecode, tstep / 60)
list_file_MetForcing = glob.glob(os.path.join(dir_start, metfile_pattern))
# load as DataFrame:
df_forcing = sp.load_SUEWS_MetForcing_df(list_file_MetForcing[1])
# load as dict (faster for simulation if performance is heavily concerned)
# dict_forcing = sp.load_SUEWS_MetForcing_dict(list_file_MetForcing[1])
# dict_forcing.keys()[-1]


# main calulation:
# compact form:
# reload(sp)
df_forcing_part = df_forcing.iloc[:8000]
dict_forcing_part = df_forcing_part.to_dict('index')
# dict_forcing_part['metforcingdata_grid'] = np.array(
#     df_forcing_part.values, dtype=np.float, order='F')
dict_state_init[1]['xwf']=0
dict_output, dict_state = sp.run_suews(
    df_forcing_part, {1: dict_state_init[1]})


# post-processing of model ouptuts:
# convert dict of raw output to easier DataFrame:
# {grid: Dataframe by group ({'SUEWS','ESTM','snow'})}
df_output = sp.pack_df_output(dict_output)
df_state = sp.pack_df_grid(dict_state)


# plot some variables
# output
xx = df_output.loc[1, 'SUEWS'].loc[:, ['QS', 'QN']].plot.line()
plt.show(xx)


# state variable
yy = df(df_state.loc[1, 'soilmoist'][:, [3, 2]]).plot()
plt.show(yy)


# load observations
list_file_Obs = glob.glob(os.path.join(
    dir_start, 'Obs/{}*txt'.format(filecode)))
res_obs = pd.read_table(list_file_Obs[0])
res_obs.rename(columns={'dectime': 'Dectime', 'QgHFP': 'QS_obs'}, inplace=True)


# comparison
# specify the variable for comparison
res_sim = df_output.loc[1, 'SUEWS'].loc[:, ['Dectime', 'QS']]
res_sim.rename(columns={'QS': 'QS_sim'}, inplace=True)
# merge two DF's
res_comp = res_obs.copy()
res_comp = res_obs.merge(res_sim, on='Dectime')
res_comp.head()
res_comp.plot(x='Dectime', y=['QS_obs', 'QS_sim'], kind='line',figsize=(10.5,3))
plt.show()




#
