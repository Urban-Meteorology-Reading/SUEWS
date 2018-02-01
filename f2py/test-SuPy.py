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
dict_init = sp.init_SUEWS_dict(dir_input)


# load met forcing
filecode = dict_init[1]['mod_config']['filecode']
resolutionfilesin = dict_init[1]['mod_config']['resolutionfilesin']
list_file_MetForcing = glob.glob(os.path.join(
    dir_input, '{}*{}*txt'.format(filecode, resolutionfilesin / 60 / 12)))
# load as DataFrame:
df_forcing = sp.load_SUEWS_MetForcing_df(list_file_MetForcing[0])
# load as dict (faster for simulation if performance is heavily concerned)
# dict_forcing = sp.load_SUEWS_MetForcing_dict(list_file_MetForcing[0])


# main calulation:
# compact form:
dict_output, dict_state = sp.run_suews(
    df_forcing.iloc[:10].T.to_dict(), dict_init)

# dict_res = sp.run_suews(dict(dict_forcing.items()[:10]), dict_init)


# line-by-line form (better control if manipulation is needed):
# # initialise dicts for holding results and model states
# dict_output = {}
# dict_state_grid = {grid: sub_dict['state_init']
#                    for grid, sub_dict in dict_init.items()}
# # temporal loop
# for tstep in df_forcing.index[:3000]:
#     # print 'cal at: ',tstep
#     # initialise output of tstep:
#     dict_output.update({tstep: {}})
#     # load met_forcing if the same across all grids:
#     met_forcing_tstep = sp.proc_met_forcing(df_forcing, tstep)
#     # spatial loop
#     for grid in dict_state_grid.keys():
#         state_old = dict_state_grid[grid]
#         mod_config = dict_init[grid]['mod_config']
#         # xx=sp.suews_cal_tstep(
#         #     state_old, met_forcing_tstep, mod_config)
#         # calculation at one step:
#         state_new, output_tstep = sp.suews_cal_tstep(
#             state_old, met_forcing_tstep, mod_config)
#         # update model state
#         dict_state_grid[grid].update(state_new)
#         # update output at tstep for the current grid
#         dict_output[tstep].update({grid: output_tstep})
#         # print 'end at ',tstep,'for',grid
#
# # convert dict_output to DataFrame
# df_res = pd.DataFrame.from_dict(dict_output)
#
# df_res_grid = pd.DataFrame(dict_output).T.stack().swaplevel()
# df_res_grid.index.get_level_values(0).unique().tolist()
# dict_grid_time = {grid: sp.pack_dict_output_grid(
#     df_res_grid[grid]) for grid in df_res_grid.index.get_level_values(0).unique().tolist()}
#
# sp.pack_dict_output_grid(df_res_grid[1])

# post-processing:
# convert dict of raw output to easier DataFrame:
# {grid: Dataframe by group ({'SUEWS','ESTM','snow'})}
df_output=sp.pack_df_output(dict_output)

xx=df_output.loc[1,'SUEWS'].loc[:,['Lup','Ldown']].plot.line()
plt.show(xx)

# sp.pack_dict_output_grid(df_res.loc[1])['SUEWS']
# # plot some variables
# len_var = len(df_res.loc[1, 'dataoutline'])
# fig, axes = plt.subplots(3, 1)
# for xfig in range(3):
#     axes[xfig].plot(df_res.loc[0:100, 'dataoutline'][xfig])
# fig
