import numpy as np
import suewsdataprocessing_v3
import suewsplotting_v1
import subprocess
import f90nml
import os
import sys

su = suewsdataprocessing_v3.SuewsDataProcessing()
pl = suewsplotting_v1.SuewsPlotting()

# read namelist, Runcontrol.nml
nml = f90nml.read('RunControl.nml')
fileinputpath = nml['runcontrol']['fileinputpath']
fileoutputpath = nml['runcontrol']['fileoutputpath']
filecode = nml['runcontrol']['filecode']
multiplemetfiles = nml['runcontrol']['multiplemetfiles']

# Working folder
wf = os.getcwd()
pf = sys.platform
if pf == 'win32':
    prog_name = 'SUEWS_V2015a.exe'
if pf == 'mac':
    prog_name = 'SUEWS_V2015a'
if pf == 'linux2':
    prog_name = 'SUEWS_V2015a'

### open SiteSelect to get year and gridnames
SiteIn = wf + fileinputpath[1:] + 'SUEWS_SiteSelect.txt'
f = open(SiteIn)
lin = f.readlines()
index = 2
loop_out = ''
YYYY_old = -9876
while loop_out != '-9':
    lines = lin[index].split()
    YYYY = int(lines[1])

    ### Create 5 min file
    if multiplemetfiles == 0: # one metfile
        if index == 2:
            gridcode1 = lines[0]
            data_in = wf + fileinputpath[1:] + filecode + gridcode1 + '_data.txt'
            met_old = np.loadtxt(data_in, skiprows=1)
            met_new = su.tofivemin_v1(met_old)
    else:  # multiple metfiles
        if index == 2:
            gridcode1 = lines[0]
            data_in = wf + fileinputpath[1:] + filecode + gridcode1 + '_data.txt'
            met_old = np.loadtxt(data_in, skiprows=1)
            met_new = su.tofivemin_v1(met_old)
        else:
            gridcode2 = lines[0]
            if gridcode2 != gridcode1:
                gridcode1 = gridcode2
                data_in = wf + fileinputpath[1:] + filecode + gridcode1 + '_data.txt'
                met_old = np.loadtxt(data_in, skiprows=1)
                met_new = su.tofivemin_v1(met_old)

    ### find start end end of 5 min file
    posstart = (np.where((met_new[:, 0] == YYYY) & (met_new[:, 1] == 1) & (met_new[:, 2] == 0) & (met_new[:, 3] == 5)))
    posend = (np.where((met_new[:, 0] == YYYY + 1) & (met_new[:, 1] == 1) & (met_new[:, 2] == 0) & (met_new[:, 3] == 0)))
    fixpos = 1
    if len(posstart[0]) == 0:
        posstart = 0
    if len(posend[0]) == 0:
        posend = met_new.shape[0]
        fixpos = 0

    met_save = met_new[posstart[0]:posend[0] + fixpos, :]

    ### save file
    data_out = wf + fileinputpath[1:] + filecode + gridcode1 + '_' + str(YYYY) + '_data_5.txt'
    header = 'iy id it imin qn qh qe qs qf U RH Tair pres rain kdown snow ldown fcld wuh xsmd lai kdiff kdir wdir'
    numformat = '%3d %2d %3d %2d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.4f %6.2f %6.2f %6.2f %6.2f ' \
                '%6.4f %6.2f %6.2f %6.2f %6.2f %6.2f'
    np.savetxt(data_out, met_save, fmt=numformat, delimiter=' ', header=header, comments='')
    f_handle = file(data_out, 'a')
    endoffile = [-9, -9]
    np.savetxt(f_handle, endoffile, fmt='%2d')
    f_handle.close()

    lines = lin[index + 1].split()
    loop_out = lines[0]
    index += 1


### This part runs the model ###
suewsstring = wf + '/' + prog_name
subprocess.call(suewsstring)


### This part makes hourly averages from SUEWS 5 min output ###

### Delete 5 min files
KeepTstepFilesIn = nml['runcontrol']['KeepTstepFilesIn']
KeepTstepFilesOut = nml['runcontrol']['KeepTstepFilesOut']

## open SUEWS_output.f95 to get format and header of output file NOT READY
# SiteIn = wf + '/SUEWS_Output.f95'
# f = open(SiteIn)
# lin = f.readlines()
# test = lin[0]
# test2 = re.match('sub', test)
TimeCol = np.array([1, 2, 3, 4, 5]) - 1
SumCol = np.array([18, 19, 20, 21, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 67, 68, 69]) - 1
LastCol = np.array([22, 23, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 64, 65, 66, 70]) - 1

header = '%iy id it imin dectime ' \
         'kdown kup ldown lup Tsurf qn h_mod e_mod qs QF QH QE ' \
         'P/i Ie/i E/i Dr/i ' \
         'St/i NWSt/i surfCh/i totCh/i ' \
         'RO/i ROsoil/i ROpipe ROpav ROveg ROwater ' \
         'AdditionalWater FlowChange WU_int WU_EveTr WU_DecTr WU_Grass ' \
         'RA RS ustar L_mod Fcld ' \
         'SoilSt smd smd_Paved smd_Bldgs smd_EveTr smd_DecTr smd_Grass smd_BSoil ' \
         'St_Paved St_Bldgs St_EveTr St_DecTr St_Grass St_BSoil St_Water ' \
         'LAI ' \
         'qn1_SF qn1_S Qm QmFreez Qmrain SWE Mw MwStore snowRem_Paved snowRem_Bldgs ChSnow/i alb_snow '

numformat = '%3d %2d %3d %2d %8.5f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'

for j in range(2, index):
    lines = lin[j].split()
    YYYY = int(lines[1])
    gridcode = lines[0]
    data_out = wf + fileinputpath[1:] + filecode + gridcode + '_' + str(YYYY) + '_data_5.txt'
    suews_5min = wf + fileoutputpath[1:] + filecode + gridcode + '_' + str(YYYY) + '_5.txt'
    suews_out = wf + fileoutputpath[1:] + filecode + gridcode + '_' + str(YYYY) + '_60.txt'
    suews_in = np.loadtxt(suews_5min, skiprows=1)
    suews_1hour = su.from5minto1hour_v1(suews_in, SumCol, LastCol, TimeCol)
    np.savetxt(suews_out, suews_1hour, fmt=numformat, delimiter=' ', header=header, comments='')  #, fmt=numformat

    if KeepTstepFilesIn == 0:
        os.remove(data_out)

    if KeepTstepFilesOut == 0:
        os.remove(suews_5min)


### plot results ###
# read namelist, plot.nml
plotnml = f90nml.read('plot.nml')
plotbasic = plotnml['plot']['plotbasic']
plotmonthlystat = plotnml['plot']['plotmonthlystat']
plotmonthlystat_col = plotnml['plot']['plotmonthlystat_col']

dectime = pl.make_dectime(suews_1hour)

if plotbasic == 1:
    pl.plot1hour(suews_1hour, met_old, dectime)

if plotmonthlystat == 1:
    pl.plotmonthlystatistics(suews_1hour, plotmonthlystat_col, dectime)


