import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib.pyplot as plt
# import collections
reload(sp)

# initialise SUEWS settings
dir_input = './Input'
dict_mod_cfg, dict_state_init = sp.init_SUEWS_dict(dir_input)

# load met forcing
filecode = dict_mod_cfg['filecode']
tstep = dict_state_init[1]['tstep']
list_file_MetForcing = glob.glob(os.path.join(
    dir_input, '{}*{}*txt'.format(filecode, tstep / 60)))
# load as DataFrame:
df_forcing = sp.load_SUEWS_MetForcing_df(list_file_MetForcing[0])
# load as dict (faster for simulation if performance is heavily concerned)
# dict_forcing = sp.load_SUEWS_MetForcing_dict(list_file_MetForcing[0])


# main calulation:
# compact form:
reload(sp)
dict_output, dict_state = sp.run_suews(
    df_forcing.iloc[:10].T.to_dict(), dict_state_init)


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
