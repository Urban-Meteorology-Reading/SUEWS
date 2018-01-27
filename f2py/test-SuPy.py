import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib.pyplot as plt
import collections
reload(sp)

# initialise SUEWS settings
dir_input = './Input'
dict_init = sp.init_SUEWS_dict(dir_input)


# load met forcing
filecode = dict_init[1]['mod_config']['filecode']
resolutionfilesin = dict_init[1]['mod_config']['resolutionfilesin']
list_file_MetForcing = glob.glob(os.path.join(
    dir_input, '{}*{}*txt'.format(filecode, resolutionfilesin / 60 / 12)))
df_forcing = sp.load_SUEWS_MetForcing_df(list_file_MetForcing[0])
# print df_forcing.columns
# df_forcing.rename(columns={
#     '%' + 'iy': 'iy',
#     'id': 'id',
#     'it': 'it',
#     'imin': 'imin',
#     'Kdn': 'avkdn',
#     'RH': 'avrh',
#     'Wind': 'avu1',
#     'fcld': 'fcld_obs',
#     'lai_hr': 'lai_obs',
#     'ldown': 'ldown_obs',
#     'rain': 'precip',
#     'press': 'press_hpa',
#     'QH': 'qh_obs',
#     'Q*': 'qn1_obs',
#     'snow': 'snow_obs',
#     'Td': 'temp_c',
#     'xsmd': 'xsmd',
#     'all': 'metforcingdata_grid'})
# print df_forcing.columns
#
#
# df_forcing[['id', 'imin']] = df_forcing[['id', 'imin']].apply(
#     lambda ss: ss.astype(int))
# df_forcing['id']


# pd.to_numeric(df_forcing['id'],downcast='int16')

df_res = sp.run_suews(df_forcing.iloc[:], dict_init)
#
# reload(sp)
# sp.proc_met_forcing(df_forcing, 100)
# sp.proc_met_forcing(df_forcing, 3333)
# # main calculation
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

# sp.pack_dict_output_grid(df_res.loc[1])['SUEWS']
# # plot some variables
# len_var = len(df_res.loc[1, 'dataoutline'])
# fig, axes = plt.subplots(3, 1)
# for xfig in range(3):
#     axes[xfig].plot(df_res.loc[0:100, 'dataoutline'][xfig])
# fig
