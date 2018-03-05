import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
# from pandas import DataFrame as df
import matplotlib.pyplot as plt
# import copy
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
# df_output = sp.pack_df_output(dict_output)
# df_state = sp.pack_df_grid(dict_state)


# %% plot some variables
idx = pd.IndexSlice
# output
xx = df_output.loc[idx[:, 1], [('datetime', 'Dectime'),
                               ('SUEWS', 'QS'), ('SUEWS', 'QN')]]
xx.columns = xx.columns.droplevel()
xx.plot(x='Dectime')
plt.show()

# state variable
idx = pd.IndexSlice

df_state_time = pd.concat(
    [df_state,
     df_output.loc[:, ('datetime', 'Dectime')].rename('Dectime')],
    axis=1).dropna()
# entry with single value
y1 = df_state_time.loc[idx[:, 1], :].plot(
    y='aerodynamicresistancemethod', x='Dectime')
plt.show()
# entry with value arrays
y2 = df_state_time.loc[:, ['Dectime', 'soilmoist']]
y2 = y2.join(y2['soilmoist'].apply(pd.Series), rsuffix='ss_').drop(
    columns=['soilmoist']).loc[:,['Dectime',1,3]].plot(x='Dectime')
plt.show()

y3 = df_state_time.loc[100:, ['Dectime', 'qn1_av_store']]
y3['qn1_av_store'] = y3['qn1_av_store'].apply(np.mean)
y3=y3.plot(x='Dectime')
plt.show()


# %% comparison
# load observations
list_file_Obs = glob.glob(
    os.path.join(dir_start,'Obs', '{}*txt'.format(filecode)))
res_obs = pd.read_table(list_file_Obs[0], delim_whitespace=True)
res_obs.columns = ['Dectime', 'QS_obs']
# res_obs.rename(columns=, inplace=True)

# specify the variable for comparison
grid = df_output.index[0][1]
res_sim = df_output.loc[idx[:, grid], [
    ('datetime', 'Dectime'), ('SUEWS', 'QS'), ('SUEWS', 'QN')]].copy()
res_sim.columns = res_sim.columns.droplevel()
res_sim.rename(columns={'QS': 'QS_sim'}, inplace=True)
# merge two DF's
res_comp = res_sim.copy()
res_comp = res_comp.merge(res_obs, left_on='Dectime',
                          right_on='Dectime', how='inner')
# res_comp.head()
plt.close('all')
res_comp.loc[:].plot(x='Dectime', y=['QS_obs', 'QS_sim'],
                     kind='line', figsize=(8, 3), style=['o'])
plt.show()

#
