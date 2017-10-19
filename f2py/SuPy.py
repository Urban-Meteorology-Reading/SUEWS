##########################################################################
# SUEWS for Python
# Ting Sun, ting.sun@reading.ac.uk
#
# Description:
# the core of this module is a f2py-based SUEWS calculation core
# with all physical functions embedded; the IO and flow control parts
# are handled by python.
#
# History:
# 18 Oct 2017: initial version
##########################################################################

# load dependency modules
import numpy as np
import pandas as pd
import os
import glob
import f90nml
# load configurations: mod_config
dir_input = './input'
os.listdir(os.path.join(dir_input))

# minimal required input files for configuration:
list_file_input = ['RunControl.nml',
                   'SUEWS_AnthropogenicHeat.txt',
                   'SUEWS_BiogenCO2.txt',
                   'SUEWS_Conductance.txt',
                   'SUEWS_ESTMCoefficients.txt',
                   'SUEWS_ESTMCoefficients_asphaltcoeffs.txt',
                   'SUEWS_Irrigation.txt',
                   'SUEWS_NonVeg.txt',
                   'SUEWS_OHMCoefficients.txt',
                   'SUEWS_Profiles.txt',
                   'SUEWS_SiteSelect.txt',
                   'SUEWS_Snow.txt',
                   'SUEWS_Soil.txt',
                   'SUEWS_Veg.txt',
                   'SUEWS_Water.txt',
                   'SUEWS_WithinGridWaterDist.txt']

# RunControl.nml
df_RunControl = pd.DataFrame(f90nml.read(os.path.join(dir_input,'RunControl.nml')))
df_RunControl.loc[:,'runcontrol'].index

df_RunControl.loc['fileinputpath','runcontrol']
for x in df_RunControl.runcontrol.iteritems():
    # print x
    if type(x[1])==str:
        print '{var}={val:{c}^{n}}'.format(var=x[0],val=x[1],n=len(x[1])+2,c='\'')
    else:
        print '{var}={val}'.format(var=x[0],val=x[1])
    # exec('{}={}'.format(*x))

sx='dog321432'
'{s:{c}^{n}}'.format(s=sx,n=len(sx)+2,c='x')


print 'x="x"'
str(x[1])

# load surface characteristics


# load meterological forcing data: met_forcing_array


# initial state: state_init


# higher-level wrapper for suews_cal_main
# [state_new,output_tstep]=suews_cal_tstep(state_old,met_forcing_tstep,mod_config)


# end
print 'good here'
