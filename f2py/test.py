#%%
import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib.pyplot as plt
import copy
# import collections
#%%
reload(sp)
# initialise SUEWS settings
# dir_input = './Input'
dir_start = '../SampleRun'
ser_mod_cfg, df_state_init = sp.init_SUEWS_pd(dir_start)
# load met forcing
filecode, kdownzen, tstep_met_in, tstep_ESTM_in, dir_input0 = ser_mod_cfg[
    ['filecode', 'kdownzen',
     'resolutionfilesin', 'resolutionfilesinestm', 'fileinputpath']]
tstep_mod, lat, lon, alt, timezone = df_state_init.iloc[0][
    ['tstep', 'lat', 'lng', 'alt', 'timezone']]
dir_input = os.path.join(dir_start, dir_input0)
grid = df_state_init.index[0]
forcingfile_met_pattern = os.path.join(
    dir_input,
    '{site}{grid}*{tstep}*txt'.format(
        site=filecode,
        grid=(grid if ser_mod_cfg['multiplemetfiles'] == 1 else ''),
        tstep=tstep_met_in / 60))
# metfile_pattern = 'Input/{}*{}*txt'.format(filecode, tstep_in / 60)
list_file_MetForcing = [
    f for f in glob.glob(forcingfile_met_pattern)
    if 'ESTM' not in f]
forcingfile_ESTM_pattern = os.path.join(
    dir_input,
    '{site}{grid}*{tstep}*txt'.format(
        site=filecode,
        grid=(grid if ser_mod_cfg['multipleestmfiles'] == 1 else ''),
        tstep=tstep_ESTM_in / 60))
list_file_ESTMForcing = [
    f for f in glob.glob(forcingfile_ESTM_pattern)
    if 'ESTM' in f]
# %%
reload(sp)
multiplemetfiles = ser_mod_cfg['multiplemetfiles']
xx_met_df = sp.load_SUEWS_Forcing_met_df_raw(
    dir_input, filecode, grid, tstep_met_in, multiplemetfiles)

multipleestmfiles = ser_mod_cfg['multipleestmfiles']
xx_estm_df = sp.load_SUEWS_Forcing_ESTM_df_raw(
    dir_input, filecode, grid, tstep_ESTM_in, multipleestmfiles)
xx_estm_df.head()
# ser_var2siteselect=pd.Series(sp.dict_var2SiteSelect)
# ser_var2siteselect.to_json('var2siteselect.json')
# import os


#%%
df_forcing = sp.load_SUEWS_MetForcing_df_resample(
    list_file_MetForcing[0],
    tstep_in, tstep_mod, lat, lon, alt, timezone, 0)

df_forcing_fort = sp.load_SUEWS_MetForcing_df_raw(
    '../SampleRun/Input/Saeve_2004_data_5.txt')
np.max(df_forcing_fort.loc[:, ['id', 'it', 'imin', 'avkdn']
                           ] - df_forcing.loc[:, ['id', 'it', 'imin', 'avkdn']])
df_forcing.loc[50:100, ['id', 'it', 'imin', 'avkdn']]

# %%
reload(sp)
grid
xx_forcing_df = sp.load_SUEWS_Forcing_df_resample(
    dir_start, grid, ser_mod_cfg, df_state_init)


xx_forcing_df.head()







#
