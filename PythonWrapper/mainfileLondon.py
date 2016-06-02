import numpy as np
import suewsdataprocessing_v3
import Tkinter
import tkFileDialog
import matplotlib.pylab as plt

su = suewsdataprocessing_v3.SuewsDataProcessing()

translate = int(input('Put data in columns (1) or do time interpolation (0)?: '))

if translate == 1:
    # This part will be used for users to get the parameters in correct columns
    old = 1 # This should be 0. Only used for old SOLWEIG-files
    ver = 2015 # Version of our model
    #Tkinter.Tk().withdraw() # Close the root window
    inputdata = tkFileDialog.askopenfilename(filetypes=[("Text Files", "*.txt;*.csv;*.prn")])
    outputdata = tkFileDialog.asksaveasfilename(filetypes=[("Text Files", "*.txt")])
    delim = ' ' # delimiter of output file
    su.translatemetdata(old, ver, inputdata, outputdata, delim)

else:
    timeinterpol = int(input('5 min to 1 hour (1) or 1 hour to 5 min (0)?: '))
    delim = ' ' # delimiter of output file
    if timeinterpol == 0: ## This part converts metdata into 5 min0
        #data_in = 'Kc1_2012_data.txt'
        #data_out = 'Kc1_2012_data_5min.txt'
        data_in = tkFileDialog.askopenfilename(filetypes=[("Text Files", "*.txt;*.csv;*.prn")])
        data_out = tkFileDialog.asksaveasfilename(filetypes=[("Text Files", "*.txt")])
        met_old = np.loadtxt(data_in, skiprows=1)

        met_new = su.tofivemin_v1(met_old)

        # Save as text files
        header = 'iy id it imin qn qh qe qs qf U RH Tair pres rain kdown snow ldown fcld wuh xsmd lai kdiff kdir wdir'
        numformat = numformat = '%3d %2d %3d %2d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.4f %6.2f %6.2f %6.2f %6.2f %6.4f %6.2f %6.2f %6.2f %6.2f %6.2f'
        np.savetxt(data_out, met_new, fmt=numformat, delimiter=delim, header=header, comments='')
        f_handle = file(data_out, 'a')
        endoffile = [-9, -9]
        np.savetxt(f_handle, endoffile, fmt='%2d')
        f_handle.close()

    else: # This part makes hourly averages from SUEWS output
        #suews_5min = 'Sm1_2010_5.txt'
        #suews_out = 'Sm1_2010_60.txt'
        suews_5min = tkFileDialog.askopenfilename(filetypes=[("Text Files", "*.txt;*.csv;*.prn")])
        suews_out = tkFileDialog.asksaveasfilename(filetypes=[("Text Files", "*.txt")])
        suews_in = np.loadtxt(suews_5min, skiprows=1)
        suews_1hour = su.from5minto1hour_v1(suews_in)

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

        # 13 15 15 15 15 15 5 15 15 15 15 15 15 8
        np.savetxt(suews_out, suews_1hour, fmt=numformat, delimiter=' ', header=header, comments='')