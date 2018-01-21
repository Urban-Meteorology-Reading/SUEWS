import SuPy.SuPy_module as sp
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline

# initialise SUEWS settings
dir_input = './Input'
dict_init = sp.init_SUEWS_dict(dir_input)


# load met forcing
filecode = dict_init[1]['mod_config']['filecode']
resolutionfilesin = dict_init[1]['mod_config']['resolutionfilesin']
list_file_MetForcing = glob.glob(os.path.join(
    dir_input, '{}*{}*txt'.format(filecode, resolutionfilesin / 60)))
df_forcing = sp.load_SUEWS_MetForcing(list_file_MetForcing[0])


# main calculation
state_old = dict_init[1]['state_init']
mod_config = dict_init[1]['mod_config']
dict_output={}
for tstep in df_forcing.index:
    met_forcing_tstep = sp.proc_met_forcing(df_forcing, 0)
    state_new, output_tstep = sp.suews_cal_tstep(
        state_old, met_forcing_tstep, mod_config)
    state_old=state_new
    dict_output.update({tstep:output_tstep})

# convert dict_output to DataFrame
df_res=pd.DataFrame.from_dict(dict_output).T

# plot some variables
len_var=len(df_res.loc[1,'dataoutline'])
fig,axes=plt.subplots(3,1);
for xfig in range(3):
    axes[xfig].plot(df_res.loc[0:100,'dataoutline'][xfig])
