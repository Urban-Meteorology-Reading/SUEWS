##############################################################################
# SuPy Demo
# SuPy: SUEWS that speaks Python
#

# Author: Ting Sun, ting.sun@reading.ac.uk
# History:
# 06 Mar 2018: alpha release

# SUEWS manual: http://urban-climate.net/umep/SUEWS
# note:
# this script demonstrates how SuPy can be used to:
# 1. perform SUEWS simulaitons
# 2. post-processing and visualisation
##############################################################################

##############################################################################
# %% load dependencies
import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# predefine index slicer for matplotlib plotting
idx = pd.IndexSlice

##############################################################################
# %% initialise SUEWS settings
# dir_start: the path to the SUEWS simulaiton direcotry where RunControl
# is placed
dir_start = '../SampleRun'

# init_SUEWS_pd initialise the SUEWS by pre-processing the input folder
# and save all configuration info into two pandas data structures:
# a: ser_mod_cfg: a Series that stores model-wide configuration info
# b: df_state_init: a DataFrame that stores grid-specific info that can be used as
#                   initial conditions for SUEWS simulations
ser_mod_cfg, df_state_init = sp.init_SUEWS_pd(dir_start)

# load met forcing as DataFrame:
# grid name to for met forcing
grid = df_state_init.index[0]
# load met forcing for the specified grid
df_forcing = sp.load_SUEWS_Forcing_df_grid(
    dir_start, grid, ser_mod_cfg, df_state_init)


# main calulation:
# here we only perform simulation for a certain number of timesteps
df_forcing_part = df_forcing.iloc[:288 * 3]
# `xwf` a test variable, leave it as it is for now: DON'T REMOVE
df_state_init.loc[:, 'xwf'] = 0
# run_suews_df performs SUEWS simulaiton and two DataFrames will be generated:
# df_output: an array of all output variables
# df_state: an array of all state variables, each of them can be used as an initial condition
df_output, df_state = sp.run_suews_df(df_forcing_part, df_state_init)


##############################################################################
# %% plot some variables
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
    columns=['soilmoist']).loc[:, ['Dectime', 1, 3]].plot(x='Dectime')
plt.show()

y3 = df_state_time.loc[100:, ['Dectime', 'qn1_av_store']]
y3['qn1_av_store'] = y3['qn1_av_store'].apply(np.mean)
y3 = y3.plot(x='Dectime')
plt.show()


# %% comparison
# load observations
filecode = ser_mod_cfg['filecode']
list_file_Obs = glob.glob(
    os.path.join(dir_start, 'Obs', '{}*txt'.format(filecode)))
res_obs = pd.read_table(list_file_Obs[0], delim_whitespace=True)
res_obs.columns = ['Dectime', 'QS_obs']

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
