##########################################################################
# SUEWS for Python
#
# Authors:
# Ting Sun, ting.sun@reading.ac.uk
# Stefan Smith, s.t.smith@reading.ac.uk
#
# Description:
# the core of this module is a f2py-based SUEWS calculation core
# with all physical functions embedded; the IO and flow control parts
# are handled by python.
#
# History:
# 18 Oct 2017, TS and SS: initial version
# 22 Oct 2017, TS: added dictionary for looking up surface Characteristics
##########################################################################

# load dependency modules
import numpy as np
import pandas as pd
import os
import glob
import f90nml
import SUEWS_driver
reload(SUEWS_driver)
from SUEWS_driver import suews_driver as sd
# load configurations: mod_config
dir_input = './input'
# os.listdir(os.path.join(dir_input))

# minimal required input files for configuration:
list_file_input = ['RunControl.nml',
                   'SUEWS_AnthropogenicHeat.txt',
                   'SUEWS_BiogenCO2.txt',
                   'SUEWS_Conductance.txt',
                   'SUEWS_ESTMCoefficients.txt',
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

# library of all properties
dict_libVar2File = {fileX.replace('.txt', '').replace(
    'SUEWS', 'lib'): fileX for fileX in list_file_input
    if fileX.endswith('.txt')}

# list of file paths
list_filepath = [os.path.join(dir_input, fileX) for fileX in list_file_input]


# this function can handle all SUEWS nml files
def load_SUEWS_nml(file):
    df = pd.DataFrame(f90nml.read(os.path.join(dir_input, 'RunControl.nml')))
    return df
    # list_cmd=[]
    # for x in df.iloc[:, 0].iteritems():
    #     # print 'set ', x[0], '=', x[1]
    #     if type(x[1]) == str:
    #         cmd = '{var}={val:{c}^{n}}'.format(
    #             var=x[0], val=x[1], n=len(x[1]) + 2, c='\'')
    #     else:
    #         cmd = '{var}={val}'.format(var=x[0], val=x[1])
    #
    #     # append cmd to final list
    #     list_cmd.append(cmd)
    #
    # return list_cmd



# RunControl.nml
lib_RunControl=load_SUEWS_nml(os.path.join(dir_input, 'RunControl.nml'))

# for cmd in lib_RunControl:
#     exec(cmd)



# load all tables (i.e., txt files)
def load_SUEWS_table(fileX):
    rawdata = pd.read_table(fileX, delim_whitespace=True,
                            comment='!', error_bad_lines=True,
                            skiprows=1, index_col=0).dropna()
    return rawdata


# load all tables into variables staring with 'lib_' and filename
for k, v in dict_libVar2File.iteritems():
    v = os.path.join(dir_input, v)
    cmd = '{var}=load_SUEWS_table({val:{c}^{n}})'.format(
        var=k, val=v, n=len(v) + 2, c='\'')
    print cmd
    exec(cmd)

# lib_BiogenCO2


############  Load Surface Characteristics from input file ###############
# Created: 19/10/2017 Stefan T. Smith
# Generic function for loading in input data and returning data to be used
# according to specified code value.
#
# dirName  - relative directory location of input file
# fileName - name of file
# headSkip - number of lines to miss before headerline
# numberOfColumns - number of columns to be used from input file
# Code - the code number used to identify which row to return. Default set
# to -10000 and if not specified all rows returned as array
def loadSurfaceCharacteristics(dirName, fileName, headSkip, numberOfColumns, Code=-10000, columnHeaderName='Code'):

    # Note the error_bad_lines = False so poor data will be dropped rather
    # than cause an exception.
    df = pd.read_csv(os.path.join(dirName, fileName), sep="\s+",
                     header=headSkip, usecols=range(numberOfColumns),
                     error_bad_lines=False)

    if (Code == -10000):
        return df.as_matrix()

    return df.loc[df[columnHeaderName] == Code].as_matrix()

############ End Load Surface Characteristics from input file #############


# dictionary for links between code in SiteSelect to properties in according tables
# this is described in SUEWS online manual:
# http://urban-climate.net/umep/SUEWS#SUEWS_SiteSelect.txt
dict_Code2File = {
    'Code_Paved': 'SUEWS_NonVeg.txt',
    'Code_Bldgs': 'SUEWS_NonVeg.txt',
    'Code_EveTr': 'SUEWS_Veg.txt',
    'Code_DecTr': 'SUEWS_Veg.txt',
    'Code_Grass': 'SUEWS_Veg.txt',
    'Code_Bsoil': 'SUEWS_NonVeg.txt',
    'Code_Water': 'SUEWS_Water.txt',
    'CondCode': 'SUEWS_Conductance.txt',
    'SnowCode': 'SUEWS_Snow.txt',
    'SnowClearingProfWD': 'SUEWS_Profiles.txt',
    'SnowClearingProfWE': 'SUEWS_Profiles.txt',
    'AnthropogenicCode': 'SUEWS_AnthropogenicHeat.txt',
    'IrrigationCode': 'SUEWS_Irrigation.txt',
    'WaterUseProfManuWD': 'SUEWS_Profiles.txt',
    'WaterUseProfManuWE': 'SUEWS_Profiles.txt',
    'WaterUseProfAutoWD': 'SUEWS_Profiles.txt',
    'WaterUseProfAutoWE': 'SUEWS_Profiles.txt',
    'WithinGridPavedCode': 'SUEWS_WithinGridWaterDist.txt',
    'WithinGridBldgsCode': 'SUEWS_WithinGridWaterDist.txt',
    'WithinGridEveTrCode': 'SUEWS_WithinGridWaterDist.txt',
    'WithinGridDecTrCode': 'SUEWS_WithinGridWaterDist.txt',
    'WithinGridGrassCode': 'SUEWS_WithinGridWaterDist.txt',
    'WithinGridUnmanBSoilCode': 'SUEWS_WithinGridWaterDist.txt',
    'WithinGridWaterCode': 'SUEWS_WithinGridWaterDist.txt',
    'Code_ESTMClass_Paved1': 'SUEWS_ESTMCoefficients.txt',
    'Code_ESTMClass_Paved2': 'SUEWS_ESTMCoefficients.txt',
    'Code_ESTMClass_Paved3': 'SUEWS_ESTMCoefficients.txt',
    'Code_ESTMClass_Bldgs1': 'SUEWS_ESTMCoefficients.txt',
    'Code_ESTMClass_Bldgs2': 'SUEWS_ESTMCoefficients.txt',
    'Code_ESTMClass_Bldgs3': 'SUEWS_ESTMCoefficients.txt',
    'Code_ESTMClass_Bldgs4': 'SUEWS_ESTMCoefficients.txt',
    'Code_ESTMClass_Bldgs5': 'SUEWS_ESTMCoefficients.txt',
    'OHMCode_SummerWet': 'SUEWS_OHMCoefficients.txt',
    'OHMCode_SummerDry': 'SUEWS_OHMCoefficients.txt',
    'OHMCode_WinterWet': 'SUEWS_OHMCoefficients.txt',
    'OHMCode_WinterDry': 'SUEWS_OHMCoefficients.txt',
    'ESTMCode': 'SUEWS_ESTMCoefficients.txt',
    'EnergyUseProfWD': 'SUEWS_Profiles.txt',
    'EnergyUseProfWE': 'SUEWS_Profiles.txt',
    'ActivityProfWD': 'SUEWS_Profiles.txt',
    'ActivityProfWE': 'SUEWS_Profiles.txt',
    'TraffProfWD': 'SUEWS_Profiles.txt',
    'TraffProfWE': 'SUEWS_Profiles.txt',
    'PopProfWD': 'SUEWS_Profiles.txt',
    'PopProfWE': 'SUEWS_Profiles.txt',
    'BiogenCO2Code': 'SUEWS_BiogenCO2.txt',
    'SoilTypeCode':'SUEWS_Soil.txt'}



# look up properties according to code in `lib` variables
def lookup_code(codeName, codeValue):
    str_lib = dict_Code2File[codeName].replace(
        '.txt', '').replace('SUEWS', 'lib')
    str_code = '{:d}'.format(int(codeValue))
    cmd = '{lib}.loc[{code}]'.format(lib=str_lib, code=str_code)
    # print cmd
    res = eval(cmd)
    return res


# get surface properties:
# if a code, look up its value in the according table;
# if a property value, retrieve it.
def get_property(varName, varVal):
    if varName in dict_Code2File.keys():
        res = lookup_code(varName, varVal)
    else:
        res = varVal
    return res


# # construct a dict with all surface properties included
# dict_SiteSelect = {ir: {k: get_property(k, v)
#                         for k, v in row.iteritems()}
#                         for ir, row in lib_SiteSelect.iterrows()}
#
# sorted(dict_SiteSelect[1].keys())
# # for k, v in lib_SiteSelect.iloc[0].iteritems():
# #     print k, v
# #     if k in dict_Code2File.keys():
# #         print k, v, lookup_code(k, v)



# variable translation as done in Fortran-SUEWS
dict_var2SiteSelect ={
     'lat': 'lat',
     'lng': 'lng',
     'timezone': 'Timezone',
     'alt': 'Alt',
     'z': 'z',
     'ah_min': {'AnthropogenicCode': ['AHMin_WD', 'AHMin_WE']},
     'ah_slope_cooling':
        {'AnthropogenicCode': ['AHSlope_Cooling_WD', 'AHSlope_Cooling_WE']},
     'ah_slope_heating':
        {'AnthropogenicCode': ['AHSlope_Heating_WD', 'AHSlope_Heating_WE']},
     'alb': [{'Code_Paved': 'AlbedoMax'},
             {'Code_Bldgs': 'AlbedoMax'},
             {'Code_EveTr': 'AlbedoMax'},
             {'Code_DecTr': 'AlbedoMax'},
             {'Code_Grass': 'AlbedoMax'},
             {'Code_Bsoil': 'AlbedoMax'},
             {'Code_Water': 'AlbedoMax'}],
     'albmax_dectr': {'Code_EveTr': 'AlbedoMax'},
     'albmax_evetr': {'Code_DecTr': 'AlbedoMax'},
     'albmax_grass': {'Code_Grass': 'AlbedoMax'},
     'albmin_dectr': {'Code_EveTr': 'AlbedoMin'},
     'albmin_evetr': {'Code_DecTr': 'AlbedoMin'},
     'albmin_grass': {'Code_Grass': 'AlbedoMin'},
     'alpha_bioco2':
         [{'Code_EveTr': {'BiogenCO2Code': 'alpha'}},
          {'Code_DecTr': {'BiogenCO2Code': 'alpha'}},
          {'Code_Grass': {'BiogenCO2Code': 'alpha'}}],
     'alpha_enh_bioco2':
        [{'Code_EveTr': {'BiogenCO2Code': 'alpha_enh'}},
         {'Code_DecTr': {'BiogenCO2Code': 'alpha_enh'}},
         {'Code_Grass': {'BiogenCO2Code': 'alpha_enh'}}],
     'alt': 'Alt',
     'baset': [{'Code_EveTr': 'BaseT'},
               {'Code_DecTr': 'BaseT'},
               {'Code_Grass': 'BaseT'}],
     'basete': [{'Code_EveTr': 'BaseTe'},
                {'Code_DecTr': 'BaseTe'},
                {'Code_Grass': 'BaseTe'}],
     'basethdd': {'AnthropogenicCode': 'BaseTHDD'},
     'beta_bioco2':
        [{'Code_EveTr': {'BiogenCO2Code': 'beta'}},
         {'Code_DecTr': {'BiogenCO2Code': 'beta'}},
         {'Code_Grass': {'BiogenCO2Code': 'beta'}}],
     'beta_enh_bioco2':
        [{'Code_EveTr': {'BiogenCO2Code': 'beta_enh'}},
         {'Code_DecTr': {'BiogenCO2Code': 'beta_enh'}},
         {'Code_Grass': {'BiogenCO2Code': 'beta_enh'}}],
     'bldgh': 'H_Bldgs',
     'capmax_dec': {'Code_DecTr': 'StorageMax'},
     'capmin_dec': {'Code_DecTr': 'StorageMin'},
     'chanohm': [{'Code_Paved': 'AnOHM_Ch'},
                 {'Code_Bldgs': 'AnOHM_Ch'},
                 {'Code_EveTr': 'AnOHM_Ch'},
                 {'Code_DecTr': 'AnOHM_Ch'},
                 {'Code_Grass': 'AnOHM_Ch'},
                 {'Code_Bsoil': 'AnOHM_Ch'},
                 {'Code_Water': 'AnOHM_Ch'}],
     'cpanohm': [{'Code_Paved': 'AnOHM_Cp'},
                 {'Code_Bldgs': 'AnOHM_Cp'},
                 {'Code_EveTr': 'AnOHM_Cp'},
                 {'Code_DecTr': 'AnOHM_Cp'},
                 {'Code_Grass': 'AnOHM_Cp'},
                 {'Code_Bsoil': 'AnOHM_Cp'},
                 {'Code_Water': 'AnOHM_Cp'}],
     'crwmax': {'SnowCode': 'CRWMax'},
     'crwmin': {'SnowCode': 'CRWMin'},
     'daywat': {'IrrigationCode':
                ['DayWat(1)',
                 'DayWat(2)',
                 'DayWat(3)',
                 'DayWat(4)',
                 'DayWat(5)',
                 'DayWat(6)',
                 'DayWat(7)']},
     'daywatper': {'IrrigationCode':
                   ['DayWatPer(1)',
                    'DayWatPer(2)',
                    'DayWatPer(3)',
                    'DayWatPer(4)',
                    'DayWatPer(5)',
                    'DayWatPer(6)',
                    'DayWatPer(7)']},
     'dectreeh': 'H_DecTr',
     'drainrt': 'LUMPS_DrRate',
     'ef_umolco2perj': {'AnthropogenicCode': 'EF_umolCO2perJ'},
     'emis': [{'Code_Paved': 'Emissivity'},
              {'Code_Bldgs': 'Emissivity'},
              {'Code_EveTr': 'Emissivity'},
              {'Code_DecTr': 'Emissivity'},
              {'Code_Grass': 'Emissivity'},
              {'Code_Bsoil': 'Emissivity'},
              {'Code_Water': 'Emissivity'}],
     'enef_v_jkm': {'AnthropogenicCode': 'EnEF_v_Jkm'},
     'evetreeh': 'H_EveTr',
     'faibldg': 'FAI_Bldgs',
     'faidectree': 'FAI_DecTr',
     'faievetree': 'FAI_EveTr',
     'faut': {'IrrigationCode': 'Faut'},
     'fcef_v_kgkm': {'AnthropogenicCode': 'FcEF_v_kgkm'},
     'frfossilfuel_heat': {'AnthropogenicCode': 'FrFossilFuel_Heat'},
     'frfossilfuel_nonheat': {'AnthropogenicCode': 'FrFossilFuel_NonHeat'},
     'g1': {'CondCode': 'G1'},
     'g2': {'CondCode': 'G2'},
     'g3': {'CondCode': 'G3'},
     'g4': {'CondCode': 'G4'},
     'g5': {'CondCode': 'G5'},
     'g6': {'CondCode': 'G6'},
     'gddfull':
         [{'Code_EveTr': 'GDDFull'},
          {'Code_DecTr': 'GDDFull'},
          {'Code_Grass': 'GDDFull'}],
     'gsmodel': {'CondCode': 'gsModel'},
     'ie_a': {'IrrigationCode': ['Ie_a1', 'Ie_a2', 'Ie_a3']},
     'ie_end': {'IrrigationCode': 'Ie_end'},
     'ie_m': {'IrrigationCode': ['Ie_m1', 'Ie_m2', 'Ie_m3']},
     'ie_start': {'IrrigationCode': 'Ie_start'},
     'internalwateruse_h': {'IrrigationCode': 'InternalWaterUse'},
     'irrfracconif': 'IrrFr_EveTr',
     'irrfracdecid': 'IrrFr_DecTr',
     'irrfracgrass': 'IrrFr_Grass',
     'kkanohm': [{'Code_Paved': 'AnOHM_Kk'},
                 {'Code_Bldgs': 'AnOHM_Kk'},
                 {'Code_EveTr': 'AnOHM_Kk'},
                 {'Code_DecTr': 'AnOHM_Kk'},
                 {'Code_Grass': 'AnOHM_Kk'},
                 {'Code_Bsoil': 'AnOHM_Kk'},
                 {'Code_Water': 'AnOHM_Kk'}],
     'kmax': {'CondCode': 'Kmax'},
     'laimax': [{'Code_EveTr': 'LAIMax'},
                {'Code_DecTr': 'LAIMax'},
                {'Code_Grass': 'LAIMax'}],
     'laimin': [{'Code_EveTr': 'LAIMin'},
                {'Code_DecTr': 'LAIMin'},
                {'Code_Grass': 'LAIMin'}],
     #  'lai_obs': '',
     'laipower': [{'Code_EveTr': ['LeafGrowthPower1', 'LeafGrowthPower2',
                                  'LeafOffPower1', 'LeafOffPower2']},
                  {'Code_DecTr': ['LeafGrowthPower1', 'LeafGrowthPower2',
                                  'LeafOffPower1', 'LeafOffPower2']},
                  {'Code_Grass': ['LeafGrowthPower1', 'LeafGrowthPower2',
                                  'LeafOffPower1', 'LeafOffPower2']}],
     'laitype': [{'Code_EveTr': 'LAIEq'},
                 {'Code_DecTr': 'LAIEq'},
                 {'Code_Grass': 'LAIEq'}],
     'lat': 'lat',
     'lng': 'lng',
     'maxconductance': [{'Code_EveTr': 'MaxConductance'},
                        {'Code_DecTr': 'MaxConductance'},
                        {'Code_Grass': 'MaxConductance'}],
     'maxqfmetab': {'AnthropogenicCode': 'MaxQFMetab'},
     'minqfmetab': {'AnthropogenicCode': 'MinQFMetab'},
     'min_res_bioco2': [{'Code_EveTr': {'BiogenCO2Code': 'min_respi'}},
                        {'Code_DecTr': {'BiogenCO2Code': 'min_respi'}},
                        {'Code_Grass': {'BiogenCO2Code': 'min_respi'}}],
     'narp_emis_snow': {'SnowCode': 'Emissivity'},
     'narp_trans_site': 'NARP_Trans',
     'ohm_coef':
         [{'Code_Paved': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
          {'Code_Bldgs': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
          {'Code_EveTr': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
          {'Code_DecTr': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
          {'Code_Grass': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
          {'Code_Bsoil': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
          {'Code_Water': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                          {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
          {'SnowCode': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
                        {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
                        {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
                        {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]}],
     'ohm_threshsw':
        [{'Code_Paved': 'OHMThresh_SW'},
         {'Code_Bldgs': 'OHMThresh_SW'},
         {'Code_EveTr': 'OHMThresh_SW'},
         {'Code_DecTr': 'OHMThresh_SW'},
         {'Code_Grass': 'OHMThresh_SW'},
         {'Code_Bsoil': 'OHMThresh_SW'},
         {'Code_Water': 'OHMThresh_SW'},
         {'SnowCode': 'OHMThresh_SW'}],
     'ohm_threshwd':
         [{'Code_Paved': 'OHMThresh_WD'},
          {'Code_Bldgs': 'OHMThresh_WD'},
          {'Code_EveTr': 'OHMThresh_WD'},
          {'Code_DecTr': 'OHMThresh_WD'},
          {'Code_Grass': 'OHMThresh_WD'},
          {'Code_Bsoil': 'OHMThresh_WD'},
          {'Code_Water': 'OHMThresh_WD'},
          {'SnowCode': 'OHMThresh_WD'}],
     'pipecapacity': 'PipeCapacity',
     'popdensdaytime': 'PopDensDay',
     'popdensnighttime': 'PopDensNight',
     'pormax_dec': {'Code_DecTr': 'PorosityMax'},
     'pormin_dec': {'Code_DecTr': 'PorosityMin'},
     'preciplimit': {'SnowCode': 'PrecipLimSnow'},
     'preciplimitalb': {'SnowCode': 'PrecipLimAlb'},
     'qf0_beu': ['QF0_BEU_WD', 'QF0_BEU_WE'],
     'qf_a': {'AnthropogenicCode': ['QF_A_WD', 'QF_A_WE']},
     'qf_b': {'AnthropogenicCode': ['QF_B_WD', 'QF_B_WE']},
     'qf_c': {'AnthropogenicCode': ['QF_C_WD', 'QF_C_WE']},
     'radmeltfact': {'SnowCode': 'RadMeltFactor'},
     'raincover': 'LUMPS_Cover',
     'rainmaxres': 'LUMPS_MaxRes',
     'resp_a': [{'Code_EveTr': {'BiogenCO2Code': 'resp_a'}},
                {'Code_DecTr': {'BiogenCO2Code': 'resp_a'}},
                {'Code_Grass': {'BiogenCO2Code': 'resp_a'}}],
     'resp_b': [{'Code_EveTr': {'BiogenCO2Code': 'resp_b'}},
                {'Code_DecTr': {'BiogenCO2Code': 'resp_b'}},
                {'Code_Grass': {'BiogenCO2Code': 'resp_b'}}],
     'runofftowater': 'RunoffToWater',
     's1': {'CondCode': 'S1'},
     's2': {'CondCode': 'S2'},
     'sathydraulicconduct':
        [{'Code_Paved': {'SoilTypeCode': 'SatHydraulicCond'}},
         {'Code_Bldgs': {'SoilTypeCode': 'SatHydraulicCond'}},
         {'Code_EveTr': {'SoilTypeCode': 'SatHydraulicCond'}},
         {'Code_DecTr': {'SoilTypeCode': 'SatHydraulicCond'}},
         {'Code_Grass': {'SoilTypeCode': 'SatHydraulicCond'}},
         {'Code_Bsoil': {'SoilTypeCode': 'SatHydraulicCond'}}],
     'sddfull': [{'Code_EveTr': 'SDDFull'},
                 {'Code_DecTr': 'SDDFull'},
                 {'Code_Grass': 'SDDFull'}],
     'sfr': ['Fr_Paved',
             'Fr_Bldgs',
             'Fr_EveTr',
             'Fr_DecTr',
             'Fr_Grass',
             'Fr_Bsoil',
             'Fr_Water'],
     'snowalbmax': {'SnowCode': 'AlbedoMax'},
     'snowalbmin': {'SnowCode': 'AlbedoMin'},
     'snowd':
         [{'Code_Paved': 'SnowLimPatch'},
          {'Code_Bldgs': 'SnowLimPatch'},
          {'Code_EveTr': 'SnowLimPatch'},
          {'Code_DecTr': 'SnowLimPatch'},
          {'Code_Grass': 'SnowLimPatch'},
          {'Code_Bsoil': 'SnowLimPatch'},
          {'Code_Bsoil': 'SnowLimPatch'}],
     'snowdensmax': {'SnowCode': 'SnowDensMax'},
     'snowdensmin': {'SnowCode': 'SnowDensMin'},
     'snowlimbuild': {'Code_Bldgs': 'SnowLimRemove'},
     'snowlimpaved': {'Code_Paved': 'SnowLimRemove'},
     'soildepth':
        [{'Code_Paved': {'SoilTypeCode': 'SoilDepth'}},
         {'Code_Bldgs': {'SoilTypeCode': 'SoilDepth'}},
         {'Code_EveTr': {'SoilTypeCode': 'SoilDepth'}},
         {'Code_DecTr': {'SoilTypeCode': 'SoilDepth'}},
         {'Code_Grass': {'SoilTypeCode': 'SoilDepth'}},
         {'Code_Bsoil': {'SoilTypeCode': 'SoilDepth'}}],
     'soilstorecap':
        [{'Code_Paved': {'SoilTypeCode': 'SoilStoreCap'}},
         {'Code_Bldgs': {'SoilTypeCode': 'SoilStoreCap'}},
         {'Code_EveTr': {'SoilTypeCode': 'SoilStoreCap'}},
         {'Code_DecTr': {'SoilTypeCode': 'SoilStoreCap'}},
         {'Code_Grass': {'SoilTypeCode': 'SoilStoreCap'}},
         {'Code_Bsoil': {'SoilTypeCode': 'SoilStoreCap'}}],
     'statelimit':
        [{'Code_Paved': 'StateLimit'},
         {'Code_Bldgs': 'StateLimit'},
         {'Code_EveTr': 'StateLimit'},
         {'Code_DecTr': 'StateLimit'},
         {'Code_Grass': 'StateLimit'},
         {'Code_Bsoil': 'StateLimit'},
         {'Code_Water': 'StateLimit'}],
     'surf':
        [[{'Code_Paved': 'StorageMin'},
          {'Code_Bldgs': 'StorageMin'},
          {'Code_EveTr': 'StorageMin'},
          {'Code_DecTr': 'StorageMin'},
          {'Code_Grass': 'StorageMin'},
          {'Code_Bsoil': 'StorageMin'},
          {'Code_Water': 'StorageMin'}],
         [{'Code_Paved': 'DrainageEq'},
          {'Code_Bldgs': 'DrainageEq'},
          {'Code_EveTr': 'DrainageEq'},
          {'Code_DecTr': 'DrainageEq'},
          {'Code_Grass': 'DrainageEq'},
          {'Code_Bsoil': 'DrainageEq'},
          {'Code_Water': 'DrainageEq'}],
         [{'Code_Paved': 'DrainageCoef1'},
          {'Code_Bldgs': 'DrainageCoef1'},
          {'Code_EveTr': 'DrainageCoef1'},
          {'Code_DecTr': 'DrainageCoef1'},
          {'Code_Grass': 'DrainageCoef1'},
          {'Code_Bsoil': 'DrainageCoef1'},
          {'Code_Water': 'DrainageCoef1'}],
         [{'Code_Paved': 'DrainageCoef2'},
          {'Code_Bldgs': 'DrainageCoef2'},
          {'Code_EveTr': 'DrainageCoef2'},
          {'Code_DecTr': 'DrainageCoef2'},
          {'Code_Grass': 'DrainageCoef2'},
          {'Code_Bsoil': 'DrainageCoef2'},
          {'Code_Water': 'DrainageCoef2'}],
         [{'Code_Paved': 'StorageMax'},
          {'Code_Bldgs': 'StorageMax'},
          {'Code_EveTr': 'StorageMax'},
          {'Code_DecTr': 'StorageMax'},
          {'Code_Grass': 'StorageMax'},
          {'Code_Bsoil': 'StorageMax'},
          {'Code_Water': 'StorageMax'}],
         [{'Code_Paved': 'StorageMin'},
          {'Code_Bldgs': 'StorageMin'},
          {'Code_EveTr': 'StorageMin'},
          {'Code_DecTr': 'StorageMin'},
          {'Code_Grass': 'StorageMin'},
          {'Code_Bsoil': 'StorageMin'},
          {'Code_Water': 'StorageMin'}]],
     'surfacearea': 'SurfaceArea',
     'tau_a': {'SnowCode': 'tau_a'},
     'tau_f': {'SnowCode': 'tau_f'},
     'tau_r': {'SnowCode': 'tau_r'},
     't_critic_cooling':
        {'AnthropogenicCode': ['TCritic_Cooling_WD', 'TCritic_Cooling_WE']},
     't_critic_heating':
        {'AnthropogenicCode': ['TCritic_Heating_WD', 'TCritic_Heating_WE']},
     'tempmeltfact': {'SnowCode': 'TempMeltFactor'},
     'th': {'CondCode': 'TL'},
     'theta_bioco2': [{'Code_EveTr': {'BiogenCO2Code': 'theta'}},
                      {'Code_DecTr': {'BiogenCO2Code': 'theta'}},
                      {'Code_Grass': {'BiogenCO2Code': 'theta'}}],
     'timezone': 'Timezone',
     'tl': {'CondCode': 'TL'},
     'trafficrate': ['TrafficRate_WD', 'TrafficRate_WE'],
     'trafficunits': {'AnthropogenicCode': 'TrafficUnits'},
     'waterdist':
        [[{'WithinGridPavedCode': 'ToPaved'},
          {'WithinGridBldgsCode': 'ToPaved'},
          {'WithinGridEveTrCode': 'ToPaved'},
          {'WithinGridDecTrCode': 'ToPaved'},
          {'WithinGridGrassCode': 'ToPaved'},
          {'WithinGridUnmanBSoilCode': 'ToPaved'}],
         [{'WithinGridPavedCode': 'ToBldgs'},
          {'WithinGridBldgsCode': 'ToBldgs'},
          {'WithinGridEveTrCode': 'ToBldgs'},
          {'WithinGridDecTrCode': 'ToBldgs'},
          {'WithinGridGrassCode': 'ToBldgs'},
          {'WithinGridUnmanBSoilCode': 'ToBldgs'}],
         [{'WithinGridPavedCode': 'ToEveTr'},
          {'WithinGridBldgsCode': 'ToEveTr'},
          {'WithinGridEveTrCode': 'ToEveTr'},
          {'WithinGridDecTrCode': 'ToEveTr'},
          {'WithinGridGrassCode': 'ToEveTr'},
          {'WithinGridUnmanBSoilCode': 'ToEveTr'}],
         [{'WithinGridPavedCode': 'ToDecTr'},
          {'WithinGridBldgsCode': 'ToDecTr'},
          {'WithinGridEveTrCode': 'ToDecTr'},
          {'WithinGridDecTrCode': 'ToDecTr'},
          {'WithinGridGrassCode': 'ToDecTr'},
          {'WithinGridUnmanBSoilCode': 'ToDecTr'}],
         [{'WithinGridPavedCode': 'ToGrass'},
          {'WithinGridBldgsCode': 'ToGrass'},
          {'WithinGridEveTrCode': 'ToGrass'},
          {'WithinGridDecTrCode': 'ToGrass'},
          {'WithinGridGrassCode': 'ToGrass'},
          {'WithinGridUnmanBSoilCode': 'ToGrass'}],
         [{'WithinGridPavedCode': 'ToBSoil'},
          {'WithinGridBldgsCode': 'ToBSoil'},
          {'WithinGridEveTrCode': 'ToBSoil'},
          {'WithinGridDecTrCode': 'ToBSoil'},
          {'WithinGridGrassCode': 'ToBSoil'},
          {'WithinGridUnmanBSoilCode': 'ToBSoil'}],
         [{'WithinGridPavedCode': 'ToWater'},
          {'WithinGridBldgsCode': 'ToWater'},
          {'WithinGridEveTrCode': 'ToWater'},
          {'WithinGridDecTrCode': 'ToWater'},
          {'WithinGridGrassCode': 'ToWater'},
          {'WithinGridUnmanBSoilCode': 'ToWater'}],
         # the last surface type is tricky: needs to determine which goes in:
         # if ToRunoff !=0, use ToRunoff, otherwise use ToSoilStore
         [{'WithinGridPavedCode': ['ToRunoff', 'ToSoilStore']},
          {'WithinGridBldgsCode': ['ToRunoff', 'ToSoilStore']},
          {'WithinGridEveTrCode': ['ToRunoff', 'ToSoilStore']},
          {'WithinGridDecTrCode': ['ToRunoff', 'ToSoilStore']},
          {'WithinGridGrassCode': ['ToRunoff', 'ToSoilStore']},
          {'WithinGridUnmanBSoilCode': ['ToRunoff', 'ToSoilStore']}]
         ],
     'wetthresh':
         [{'Code_Paved': 'WetThreshold'},
          {'Code_Bldgs': 'WetThreshold'},
          {'Code_EveTr': 'WetThreshold'},
          {'Code_DecTr': 'WetThreshold'},
          {'Code_Grass': 'WetThreshold'},
          {'Code_Bsoil': 'WetThreshold'},
          {'Code_Water': 'WetThreshold'}],
     'year': 'Year',
     'z': 'z'}


# expand dict_Code2File for retrieving surface characteristics
dict_varSiteSelect2File = {
    x: 'SUEWS_SiteSelect.txt' for x in dict_var2SiteSelect.keys()}
dict_Code2File.update(dict_varSiteSelect2File)


# look up properties according to code
def lookup_code_sub(codeName, codeKey, codeValue):
    str_lib = dict_Code2File[codeName].replace(
        '.txt', '').replace('SUEWS', 'lib')
    str_code = '{:d}'.format(int(codeValue))
    cmd = '{lib}.loc[{code},{key:{c}^{n}}]'.format(
        lib=str_lib, code=str_code, key=codeKey, n=len(codeKey) + 2, c='\'')
    # print cmd
    res = eval(cmd)
    return res


# a recursive function to retrieve value based on key sequences
def lookup_KeySeq(indexKey, subKey, indexCode):
    # print indexKey, subKey, indexCode
    if type(subKey) is str:
        res = lookup_code_sub(indexKey, subKey, indexCode)
    elif type(subKey) is dict:
        indexKeyX, subKeyX = subKey.items()[0]
        indexCodeX = lookup_code_sub(indexKey, indexKeyX, indexCode)
        res = lookup_KeySeq(indexKeyX, subKeyX, indexCodeX)
    elif type(subKey) is list:
        res = []
        for subKeyX in subKey:
            indexCodeX = indexCode
            resX = lookup_KeySeq(indexKey, subKeyX, indexCodeX)
            res.append(resX)
    # final result
    return res


# construct a dictionary in the form: {grid:{var:value,...}}
dict_gridSurfaceChar={
    # grid:{k:v
    grid: {k:lookup_KeySeq(k,v,grid)
    for k,v in dict_var2SiteSelect.iteritems() }
    for grid in lib_SiteSelect.index}


# convert the above dict to DataFrame
df_gridSurfaceChar=pd.DataFrame.from_dict(dict_gridSurfaceChar).T

df_gridSurfaceChar.shape



# get the information of input variables for SUEWS_driver
docLines=np.array(
    sd.suews_cal_main.__doc__.splitlines(),
    dtype=str)
pos=np.where(
    np.logical_or(
    docLines=='Parameters',docLines=='Other Parameters'))
varInputLines=docLines[pos[0][0]+2:pos[0][1]-1]
varInputInfo=np.array([[xx.rstrip() for xx in x.split(':')] for x in varInputLines])
len(varInputInfo)

# variables to be fed in:
list_varToDo=np.array(list(set(varInputInfo[:,0])-set(df_gridSurfaceChar.columns)))

# method-related options:
list_varMethod=[x for x in varToDo if 'method' in x.lower()]

# TODO: profile-related:
list_varTstep=[x for x in varToDo if '_tstep' in x.lower()]

# time-related:
list_varTime=[
    'iy','dayofweek','imin','id','nsh_real','halftimestep',
    'dls','it','tstep_real','dectime']

# met-forcing-related:
list_varMetForcing=[
    'ldown_obs','avrh','press_hpa','lai_obs','precip',
    'avu1','ts5mindata_ir','avkdn','qh_obs','qn1_obs',
    'temp_c','metforcingdata']

# DailyState-related:
list_varDailyState=[
    'albdectr','albevetr',,'albgrass','gdd','hdd']

# setting-related:
list_varSetting=[
    'laicalcyes','numcapita','ity','diagqs','diagqn']

# model-tstep-state-related:
list_varModelState=[
    'qn1_av_store','icefrac','tsurf_ind','qn1_s_av_store',
    'qn1_store','tair24hr']

# output-related:
list_varModelState=[
    'dataoutestm','dataout','tsurf_ind','qn1_s_av_store']





# load meterological forcing data: met_forcing_array
# filecode gives the name to locate met forcing file
filecode

resolutionfilesin
list_file_MetForcing = glob.glob(os.path.join(
    'input', '*{}*{}*txt'.format(filecode, resolutionfilesin / 60)))


def func_parse_date(year, doy, hour, sec):
    dt = datetime.datetime.strptime(
        ' '.join([year, doy, hour, sec]), '%Y %j %H %M')
    # dt_base = datetime.datetime(int(year), 1, 1)
    # dt_delta = datetime.timedelta(int(doy),
    #                               np.dot(np.array([hour, sec],
    #                                   dtype=np.float),
    #                               [3600, 1]))
    # dt = dt_base + dt_delta
    return dt


def load_SUEWS_MetForcing(fileX):
    rawdata = pd.read_table(fileX, delim_whitespace=True,
                            comment='!', error_bad_lines=True,
                            parse_dates={'datetime': [0, 1, 2, 3]}, keep_date_col=True,
                            date_parser=func_parse_date).dropna()
    return rawdata


xx = load_SUEWS_MetForcing(list_file_MetForcing[0])
for x in xx.iloc[1, 1:]:
    print type(x)
# initial state: state_init


# higher-level wrapper for suews_cal_main
# [state_new,output_tstep]=suews_cal_tstep(state_old,met_forcing_tstep,mod_config)


# end
print 'good here'
