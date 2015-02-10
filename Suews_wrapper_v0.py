import numpy as np
import seuwsdataprocessing
import subprocess
su = seuwsdataprocessing.SuewsDataProcessing()

# Working folder
wf = 'M:/SUEWS/Vs87/'
progname = 'SUEWS_V2013a.exe'

### This part converts metdata into 5 min ###
data_in = wf + 'Input/Kc1_2012_data.txt'
data_out = wf + 'Input/Kc1_2012_data_5min.txt'
met_old = np.loadtxt(data_in, skiprows=1)
met_new = su.tofivemin_v1(met_old)
header = 'iy id it imin qn qh qe qs qf U RH Tair pres rain kdown snow ldown fcld wuh xsmd lai kdiff kdir wdir'
numformat = '%3d %2d %3d %2d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f ' \
            '%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f'
np.savetxt(data_out, met_new, fmt=numformat, delimiter=' ', header=header, comments='')
f_handle = file(data_out, 'a')
endoffile = [-9, -9]
np.savetxt(f_handle, endoffile, fmt='%2d')
f_handle.close()


### This part runs the model ###
suewsstring = wf + progname
subprocess.call(suewsstring)


### This part makes hourly averages from SUEWS 5 min output ###
suews_5min = wf + 'Output/Sm1_2010_5.txt'
suews_out = wf + 'Output/Sm1_2010_60.txt'
suews_in = np.loadtxt(suews_5min, skiprows=1)
suews_1hour = su.from5minto1hour_v1(suews_in)

header = '%iy  id   it imin dectime    kdown     kup    ldown     lup     Tsurf     qn      h_mod    e_mod     ' \
         'qs        QF       QH       QE      P/i      Ie/i     E/i      DR/i    Ch/i     ST/i    ROsoil/i  RO/i    ' \
         'ROpipe   ROpav     ROveg   ROwater    RA      RS     ustar    L_mod SoilSt_pav SoilSt_bldg SoilSt_ET ' \
         'SoilSt_DT SoilSt_IG SoilSt_UG St_pav St_bldg   St_ET    St_DT    St_IG    St_UG   St_wtr     Fcld ' \
         'SoilState    smd       LAI       Fw     addWater Iegrass/i Ietrees/i qn1_SF    qn1_S     Qm    delta_QS    ' \
         'Qrain     SWE    MwStore snowRem_pav snowRem_bldg ChSnow/i kup_pav   kup_blgskup_ET    kup_dT    kup_IG    ' \
         'kup_UG   kup_wtr   lup_pav  lup_bldg    lup_ET    lup_DT   lup_IG    lup_UG    lup_wtr    Ts_pav    ' \
         'Ts_bldg    Ts_ET     Ts_DT     Ts_IG     Ts_UG    Ts_wtr    qn_pav   qn_bldg     qn_ET     qn_DT     ' \
         'qn_IG    qn_UG     qn_wtrSWE_pav SWE_bldg   SWE_ET   SWE_DT   SWE_IG   SWE_UG  SWE_wtr SnowRem_pav SnowRem_' \
         'bldg Mw  Mw_pav  Mw_bldg    Mw_ET    Mw_DT    Mw_IG    Mw_UG   Mw_wtr     Qm   Qm_pav Qm_bldg   Qm_ET   ' \
         'Qm_DT   Qm_IG   Qm_UG  Qm_wtr  Qa_pav Qa_bldg   Qa_ET   Qa_DT   Qa_IG   Qa_UG  Qa_wtr QmFr_pav QmFr_bldg ' \
         'QmFr_ET QmFr_DT QmFr_IG QmFr_UG QmFr_wtr fr_pav fr_bldg fr_ET  fr_DT   fr_IG   fr_UG alb_snow RainSn_pav ' \
         'RainSn_bldg RainSn_ET RainSn_DT RainSn_IG RainSn_UG RainSn_wtr Qn_pavSnow Qn_blgsSnow Qn_ETSnow Qn_DTSnow ' \
         'Qn_IGSnow Qn_UG Qs_wtrSnow kup_pavSnow kup_blgsSnow kup_ETSnow kup_DTSnow kup_IGSnow kup_UGSnow kup_' \
         'wtrSnow frMelt_pav frMelt_bldg frMelt_ET frMelt_DT frMelt_IG frMelt_UG frMelt_wtrMwStore_pav MwStore_bldg ' \
         'MwStore_ET MwStore_DT MwStore_IG MwStore_UG MwStore_wtrdensSnow_pav densSnow_bldg densSnow_ET densSnow_DT ' \
         'densSnow_IG densSnow_UG densSnow_wtrSd_pav Sd_bldg Sd_ET Sd_DT Sd_IG Sd_UG Sd_water Tsnow_pav Tsnow_bldg ' \
         'Tsnow_ET Tsnow_DT Tsnow_IG Tsnow_UG Tsnow_wtr'
numformat = '%3d %2d %3d %2d %8.5f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %8.2f %8.2f %8.2f %8.2f %8.2f' \
            ' %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f' \
            ' %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f' \
            ' %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f' \
            ' %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f' \
            ' %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f' \
            ' %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f' \
            ' %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f'

# 13 15 15 15 15 15 5 15 15 15 15 15 15 8
np.savetxt(suews_out, suews_1hour, fmt=numformat, delimiter=' ', header=header, comments='')#, fmt=numformat








