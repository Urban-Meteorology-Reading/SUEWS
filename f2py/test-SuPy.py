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
ser_mod_cfg, df_state_init = sp.init_SUEWS_pd(dir_start)
# load met forcing
filecode, kdownzen, tstep_in = ser_mod_cfg[
    ['filecode', 'kdownzen', 'resolutionfilesin']]
tstep_mod, lat, lon, alt, timezone = df_state_init.iloc[0][
    ['tstep', 'lat', 'lng', 'alt', 'timezone']]
metfile_pattern = 'Input/{}*{}*txt'.format(filecode, tstep_in / 60)
list_file_MetForcing = [f
                        for f in glob.glob(os.path.join(
                            dir_start, metfile_pattern))
                        if 'ESTM' not in f]
# load as DataFrame:
df_forcing = sp.load_SUEWS_MetForcing_df_resample(
    list_file_MetForcing[0],
    tstep_in, tstep_mod, lat, lon, alt, timezone,kdownzen)
# load as dict (faster for simulation if performance is heavily concerned)
# dict_forcing = sp.load_SUEWS_MetForcing_dict(list_file_MetForcing[1])
# dict_forcing.keys()[-1]


# main calulation:
# compact form:
# reload(sp)
df_forcing_part = df_forcing.iloc[:288 * 1]
df_state_init.loc[:, 'xwf'] = 0
df_output, df_state = sp.run_suews_df(df_forcing_part, df_state_init)


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
res_comp.plot(x='Dectime', y=['QS_obs', 'QS_sim'],
              kind='line', figsize=(10.5, 3))
plt.show()


#
