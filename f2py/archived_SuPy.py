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
# 12 Jan 2018, TS: a working version for OHM with arguments sorted out
##########################################################################

# load dependency modules
import os
# dir_path = os.path.dirname(os.path.realpath(__file__))
# # os.chdir(os.getcwd()+'/f2py')
# os.chdir(dir_path)
import numpy as np
import pandas as pd
import glob
import f90nml
import datetime
# import SUEWS_driver
# reload(SUEWS_driver)
from SUEWS_driver import suews_driver as sd
from scipy import interpolate


def insensitive_glob(pattern):
    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(either, pattern)))


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
def load_SUEWS_nml(xfile):
    df = pd.DataFrame(f90nml.read(xfile))
    return df


# RunControl.nml
lib_RunControl = load_SUEWS_nml(os.path.join(dir_input, 'RunControl.nml'))

for var in lib_RunControl.index:
    val = lib_RunControl.loc[var, 'runcontrol']
    if type(val) == str:
        cmd = '{var}={val:{c}^{n}}'.format(
            var=var, val=val, n=len(val) + 2, c='\'')
    else:
        cmd = '{var}={val}'.format(var=var, val=val)
    # print cmd
    exec(cmd)


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
    # print cmd
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
    'WaterUseProfManuWD': 'SUEWS_Profiles.txt',
    'WaterUseProfManuWE': 'SUEWS_Profiles.txt',
    'WaterUseProfAutoWD': 'SUEWS_Profiles.txt',
    'WaterUseProfAutoWE': 'SUEWS_Profiles.txt',
    'BiogenCO2Code': 'SUEWS_BiogenCO2.txt',
    'SoilTypeCode': 'SUEWS_Soil.txt'}


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



# variable translation as done in Fortran-SUEWS
dict_var2SiteSelect = {
    'lat': 'lat',
    'lng': 'lng',
    'timezone': 'Timezone',
    'alt': 'Alt',
    'z': 'z',
    'ahprof':
        {'AnthropogenicCode': [
            {'EnergyUseProfWD': ':'}, {'EnergyUseProfWE': ':'}]},
    'popprof':
        {'AnthropogenicCode': [{'PopProfWD': ':'}, {'PopProfWE': ':'}]},
    'traffprof':
        {'AnthropogenicCode': [{'TraffProfWD': ':'}, {'TraffProfWE': ':'}]},
    'humactivity':
        {'AnthropogenicCode': [
            {'ActivityProfWD': ':'}, {'ActivityProfWE': ':'}]},
    'wuprofa':
        [{'WaterUseProfAutoWD': ':'}, {'WaterUseProfAutoWE': ':'}],
    'wuprofm':
        [{'WaterUseProfManuWD': ':'}, {'WaterUseProfManuWE': ':'}],
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
    'flowchange': 'FlowChange',
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
            {'Code_Bsoil': {'SoilTypeCode': 'SatHydraulicCond'}},
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
            {'Code_Bsoil': {'SoilTypeCode': 'SoilDepth'}},
            {'Code_Bsoil': {'SoilTypeCode': 'SoilDepth'}}],
    'soilstorecap':
    [{'Code_Paved': {'SoilTypeCode': 'SoilStoreCap'}},
            {'Code_Bldgs': {'SoilTypeCode': 'SoilStoreCap'}},
            {'Code_EveTr': {'SoilTypeCode': 'SoilStoreCap'}},
            {'Code_DecTr': {'SoilTypeCode': 'SoilStoreCap'}},
            {'Code_Grass': {'SoilTypeCode': 'SoilStoreCap'}},
            {'Code_Bsoil': {'SoilTypeCode': 'SoilStoreCap'}},
            {'Code_Bsoil': {'SoilTypeCode': 'SoilStoreCap'}}],
    'startDLS': 'StartDLS',
    'endDLS': 'EndDLS',
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
    if codeKey == ':':
        cmd = '{lib}.loc[{code},:].tolist()'.format(lib=str_lib, code=str_code)
    else:
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
dict_gridSurfaceChar = {
    # grid:{k:v
    grid: {k: lookup_KeySeq(k, v, grid)
           for k, v in dict_var2SiteSelect.iteritems()}
    for grid in lib_SiteSelect.index}


# convert the above dict to DataFrame
df_gridSurfaceChar = pd.DataFrame.from_dict(dict_gridSurfaceChar).T

# df_gridSurfaceChar.shape


# get the information of input variables for SUEWS_driver
docLines = np.array(
    sd.suews_cal_main.__doc__.splitlines(),
    dtype=str)
posInput = np.where(
    np.logical_or(
        docLines == 'Parameters', docLines == 'Returns'))
varInputLines = docLines[posInput[0][0] + 2:posInput[0][1] - 1]
varInputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                         for x in varInputLines])
dict_InputInfo = {xx[0]: xx[1] for xx in varInputInfo}
# print len(dict_InputInfo)

# variables to be fed in:
list_varToDo = np.array(
    list(set(varInputInfo[:, 0]) - set(df_gridSurfaceChar.columns)))

# method-related options:
list_varMethod = [x for x in list_varToDo if 'method' in x.lower()]

# profile-related:
list_varTstep = [x for x in list_varToDo if '_tstep' in x.lower()]

# time-related:
list_varTime = [
    'iy', 'dayofweek', 'imin', 'id', 'nsh_real', 'halftimestep',
    'dls', 'it', 'tstep_real', 'dectime']

# met-forcing-related:
list_varMetForcing = [
    'ldown_obs', 'avrh', 'press_hpa', 'lai_obs', 'precip',
    'avu1', 'ts5mindata_ir', 'avkdn', 'qh_obs', 'qn1_obs',
    'temp_c', 'metforcingdata', 'snow_obs', 'fcld_obs', 'xsmd']

# DailyState-related:
list_varDailyState = [
    'albdectr', 'albevetr', 'albgrass', 'gdd', 'hdd', 'lai',
    'decidcap', 'porosity', 'snowalb', 'snowdens', ]

# setting-related:
list_varSetting = [
    'laicalcyes', 'numcapita', 'ity', 'diagqs', 'diagqn']

# model-tstep-state-related:
list_varModelState = [
    'qn1_av_store', 'icefrac', 'qn1_s_av_store',
    'qn1_store', 'tair24hr', 'qn1_s_store', 'meltwaterstore',
    'snowfrac', 'snowpack', 'soilmoist', 'state']

# output-related:
list_varModelOut = [
    'dataoutestm', 'dataout']

# model-control-related:
list_varModelCtrl = [
    'ir', 'gridiv', 'veg_type',
    'ncolumnsdataout', 'readlinesmetdata', 'numberofgrids']

# ncolumnsdataout = 84  # set according to SUEWS V2017c

# extra surface characteristics variables in addition to df_gridSurfaceChar:
list_varGridSurfaceCharX = [
    'nonwaterfraction', 'pervfraction', 'vegfraction']
# a complete list of surface characteristics variables
list_varGridSurfaceChar = (
    df_gridSurfaceChar.columns.tolist()
    + list_varGridSurfaceCharX)

# load values into list_varGridSurfaceChar
list_grid = df_gridSurfaceChar.index  # list of grid
numberofgrids = len(list_grid)
for grid in list_grid:  # this is only for testing;
    # need to be moved to main loop
    for var in df_gridSurfaceChar.columns.tolist():
        # print var
        val = df_gridSurfaceChar.loc[grid, var]
        # print type(val)
        if type(val) == str:
            # print val
            cmd = '{var}={val:{c}^{n}}'.format(
                var=var, val=val, n=len(val) + 2, c='\'')
        else:
            # valX=list(val)
            # print valX
            cmd = '{var}={val}'.format(var=var, val=val)
        # print cmd
        exec(cmd)

# modify some variables by correcting dimensions
# transpoe laipower:
laipower = np.array(laipower).T
# select non-zero values for waterdist of water surface:
waterdist_water = np.array(waterdist[-1])
waterdist_water = waterdist_water[np.nonzero(waterdist_water)]
waterdist[-1] = waterdist_water
# surf order as F:
surf = np.array(surf, order='F')
# convert to np.array
alb = np.array(alb)

# derive surface fractions
# PavSurf = 1 - 1
# BldgSurf = 2 - 1
# ConifSurf = 3 - 1
DecidSurf = 4 - 1
# GrassSurf = 5 - 1
# BSoilSurf = 6 - 1
# WaterSurf = 7 - 1
# vegfraction = (sfr[ConifSurf] + sfr[DecidSurf] + sfr[GrassSurf])
# impervfraction = sfr[PavSurf] + sfr[BldgSurf]
# pervfraction = (sfr[ConifSurf] + sfr[DecidSurf]
#                 + sfr[GrassSurf] + sfr[BSoilSurf] + sfr[WaterSurf])
# nonwaterfraction = 1 - sfr[WaterSurf]

# profiles:
t_tstep = np.linspace(0, 24, num=3600 / tstep * 24, endpoint=False)
for var in list_varTstep:
    var0 = var.replace('_tstep', '')
    var0 = eval('np.array({var0}).T'.format(var0=var0))
    var0 = np.vstack((var0, var0))
    # interpolator:
    f = interpolate.interp1d(np.arange(0, 48), var0, axis=0)
    cmd = '{var}=f(t_tstep)'.format(var=var)
    # print cmd
    exec(cmd)


# now we proceed to construct several arrays for suews_cal_main
# mod_state: runtime states as initial conditions for successive steps
list_var_mod_state = (
    list_varModelState
    + list_varDailyState
    + list_varModelOut)  # list_varModelOut is in fact used as model states


# mod_forcing: forcing condition and other control variables needed for
# each time step
list_var_mod_forcing = list_varMetForcing


# mod_config: configuration info for each run/simulation
list_var_mod_config = (
    list_varModelCtrl
    + list(lib_RunControl.index)
    + list_varGridSurfaceChar)


# mod_output: outputs of each time step
# this should be parsed from suews_cal_main.__doc__
# get the information of output variables for SUEWS_driver
docLines = np.array(
    sd.suews_cal_main.__doc__.splitlines(),
    dtype=str)
posOutput = np.where(docLines == 'Returns')
varOutputLines = docLines[posOutput[0][0] + 2:]
varOutputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                          for x in varOutputLines])
dict_OutputInfo = {xx[0]: xx[1] for xx in varOutputInfo}
# len(varOutputInfo)

# transfer variable names from varOutputInfo to list_varModelOut
list_var_mod_output = varOutputInfo[:, 0]


# load meterological forcing data: met_forcing_array
# filecode and resolutionfilesin gives the name to locate met forcing file
# filecode
# resolutionfilesin
list_file_MetForcing = glob.glob(os.path.join(
    'Input', '{}*{}*txt'.format(filecode, resolutionfilesin / 60)))


def func_parse_date(year, doy, hour, sec):
    dt = datetime.datetime.strptime(
        ' '.join([year, doy, hour, sec]), '%Y %j %H %M')
    return dt


def load_SUEWS_MetForcing(fileX):
    rawdata = pd.read_table(fileX, delim_whitespace=True,
                            comment='!', error_bad_lines=True,
                            parse_dates={'datetime': [0, 1, 2, 3]},
                            keep_date_col=True,
                            date_parser=func_parse_date).dropna()
    return rawdata


df_forcing = load_SUEWS_MetForcing(list_file_MetForcing[0])
# for x in df_forcing.iloc[1, 1:]:
#     print type(x)
# initial state: state_init


# higher-level wrapper for suews_cal_main
# [state_new,output_tstep]=suews_cal_tstep(state_old,met_forcing_tstep,mod_config)

# some constant values
ndays = 366
aerodynamicresistancemethod = 2
# halftimestep = tstep / 2.
nsh = 3600 / tstep
# nsh_real = nsh * 1.
ity = 2
laicalcyes = 1
tstep_real = tstep * 1.
veg_type = 1
diagqn = 0
diagqs = 0

# initialise model state
gdd = 0 * np.ones((ndays + 1, 5), order='F')
hdd = 0 * np.ones((ndays + 5, 6), order='F')
icefrac = 0 * np.ones(7)
lai = 0. * np.ones((ndays + 5, 3), order='F')
meltwaterstore = 0 * np.ones(7, order='F')
numcapita = (popdensdaytime + popdensnighttime) / 2.
porosity = pormax_dec * np.ones(ndays + 1, order='F')
qn1_av_store = -999. * np.ones(2 * nsh + 1)
qn1_s_av_store = -999. * np.ones(2 * nsh + 1)
qn1_s_store = -999. * np.ones(nsh)
qn1_store = -999. * np.ones(nsh)
snowalb = 0. * np.ones(7)
snowdens = 0. * np.ones(7)
snowfrac = 0. * np.ones(7)
snowpack = 0. * np.ones(7)
wu_day = 0. * np.ones(((ndays + 1), 9), order='F')
soilmoist = 0.5 * np.ones(7)
state = -999. * np.ones(7)
tair24hr = 273.15 * np.ones(24 * nsh)
dayofweek = np.ones((ndays + 1, 3), order='F', dtype=np.int32)
albdectr = albmax_dectr * np.ones(ndays + 1, order='F')
albevetr = albmax_evetr * np.ones(ndays + 1, order='F')
albgrass = albmax_grass * np.ones(ndays + 1, order='F')
decidcap = surf[DecidSurf][5 - 1] * np.ones(ndays + 1, order='F')



# begin running:
nlines = len(df_forcing)
# loop over time
for ir in df_forcing.index:
    # ir = 1

    # time-related:
    dectime = (float(df_forcing.loc[ir, 'id'])
               + float(df_forcing.loc[ir, 'it']) / 24.
               + float(df_forcing.loc[ir, 'imin']) / 60.)

    # print dectime
    iy = int(df_forcing.loc[ir, '%' + 'iy'])
    id = int(df_forcing.loc[ir, 'id'])
    it = int(df_forcing.loc[ir, 'it'])
    imin = int(df_forcing.loc[ir, 'imin'])
    id_prev_t = id
    iy_prev_t = iy
    # dayofweek:
    dayofweek[id, 0] = df_forcing.loc[ir, 'datetime'].dayofweek
    dayofweek[id, 1] = df_forcing.loc[ir, 'datetime'].month
    # season: 1 for summer; 0 for winter.
    dayofweek[id, 2] = (1 if 3 < dayofweek[id, 1] < 10 else 0)
    # specify dls:
    dls = (1 if startDLS < id < endDLS else 0)

    # met forcing variables:
    metforcingdata_grid = np.array(df_forcing.iloc[:,1:],dtype=np.float)
    ts5mindata_ir = df_forcing.loc[ir, 'Td']
    avkdn = df_forcing.loc[ir, 'Kdn']
    avrh = df_forcing.loc[ir, 'RH']
    avu1 = df_forcing.loc[ir, 'Wind']
    fcld_obs = df_forcing.loc[ir, 'fcld']
    lai_obs = df_forcing.loc[ir, 'lai_hr']
    ldown_obs = df_forcing.loc[ir, 'ldown']
    precip = df_forcing.loc[ir, 'rain']
    press_hpa = df_forcing.loc[ir, 'press'] * 10
    qh_obs = df_forcing.loc[ir, 'QH']
    qn1_obs = df_forcing.loc[ir, 'Q*']
    snow_obs = df_forcing.loc[ir, 'snow']
    temp_c = df_forcing.loc[ir, 'Td']
    xsmd = df_forcing.loc[ir, 'xsmd']

    # loop over grid
    for gridiv in df_gridSurfaceChar.index:
        # gridiv = 1

        # main calculation:
        datetimeline, dataoutline, dataoutlinesnow, dataoutlineestm = \
            sd.suews_cal_main(aerodynamicresistancemethod, ah_min,
            ahprof_tstep, ah_slope_cooling, ah_slope_heating,
            alb, albdectr, albevetr, albgrass,
            albmax_dectr, albmax_evetr, albmax_grass,
            albmin_dectr, albmin_evetr, albmin_grass,
            alpha_bioco2, alpha_enh_bioco2, alt, avkdn, avrh, avu1,
            baset, basete, basethdd, beta_bioco2, beta_enh_bioco2, bldgh,
            capmax_dec, capmin_dec, chanohm, cpanohm, crwmax, crwmin,
            dayofweek, daywat, daywatper, decidcap, dectime, dectreeh,
            diagnose, diagqn, diagqs, dls, drainrt, ef_umolco2perj,
            emis, emissionsmethod, enef_v_jkm, evetreeh, faibldg,
            faidectree, faievetree, faut, fcef_v_kgkm, fcld_obs, flowchange,
            frfossilfuel_heat, frfossilfuel_nonheat, g1, g2, g3, g4, g5, g6,
            gdd, gddfull, gridiv, gsmodel, hdd, humactivity_tstep, icefrac,
            id, id_prev_t, ie_a, ie_end, ie_m, ie_start, imin,
            internalwateruse_h, irrfracconif, irrfracdecid, irrfracgrass,
            it, ity, iy, iy_prev_t, kkanohm, kmax, lai, laicalcyes,
            laimax, laimin, lai_obs, laipower, laitype, lat, ldown_obs,
            lng, maxconductance, maxqfmetab, meltwaterstore,
            metforcingdata_grid, minqfmetab,min_res_bioco2,
            narp_emis_snow, narp_trans_site, netradiationmethod, numcapita,
            ohm_coef, ohmincqf, ohm_threshsw, ohm_threshwd, pipecapacity,
            popdensdaytime, popdensnighttime, popprof_tstep,
            pormax_dec, pormin_dec, porosity, precip,
            preciplimit, preciplimitalb, press_hpa, qf0_beu,
            qf_a, qf_b, qf_c, qh_obs, qn1_av_store, qn1_obs,
            qn1_s_av_store, qn1_s_store, qn1_store, radmeltfact,
            raincover, rainmaxres, resp_a, resp_b,
            roughlenheatmethod, roughlenmommethod, runofftowater, s1, s2,
            sathydraulicconduct, sddfull, sfr, smdmethod,
            snowalb, snowalbmax, snowalbmin, snowd,
            snowdens, snowdensmax, snowdensmin, snowfrac,
            snowlimbuild, snowlimpaved, snow_obs, snowpack, snowuse,
            soildepth, soilmoist, soilstorecap, stabilitymethod,
            state, statelimit, storageheatmethod, surf, surfacearea,
            tair24hr, tau_a, tau_f, tau_r, t_critic_cooling, t_critic_heating,
            temp_c, tempmeltfact, th, theta_bioco2, timezone, tl,
            trafficrate, trafficunits, traffprof_tstep, ts5mindata_ir,
            tstep, veg_type, waterdist, waterusemethod, wetthresh,
            wu_day, wuprofa_tstep, wuprofm_tstep, xsmd, year, z)

        # print dectime, dataoutline[4:10]

        # end grid loop
    # end time loop
# end



##############################################################################
# DEPRECATED CODE
# some `historic` stuff used in the development and may be removed soemtime
# # pack up output of one grid of all tsteps
# def pack_dict_output_grid(df_grid):
#     # merge dicts of all tsteps
#     dict_group = collections.defaultdict(list)
#     for d in df_grid:
#         for k, v in d.iteritems():  # d.items() in Python 3+
#             dict_group[k].append(v)
#     # pick groups except for `datetimeline`
#     group_out = (group for group in dict_group.keys()
#                  if not group == 'datetimeline')
#     # initialise dict for holding packed output
#     dict_output_group = {}
#     # group names in lower case
#     var_df_lower = {group.lower(): group for group in var_df.index.levels[0]}
#     # pack up output of all tsteps into output groups
#     for group_x in group_out:
#         # get correct group name by cleaning and swapping case
#         group = group_x.replace('dataoutline', '').replace('line', '')
#         # print group
#         group = var_df_lower[group]
#         header_group = np.apply_along_axis(
#             list, 0, var_df.loc[['datetime', group]].index.values)[:, 1]
#         # print 'header_group', header_group
#         df_group = pd.DataFrame(
#             np.hstack((dict_group['datetimeline'], dict_group[group_x])),
#             columns=header_group)
#         # df_group[[
#         #     'Year', 'DOY', 'Hour', 'Min']] = df_group[[
#         #         'Year', 'DOY', 'Hour', 'Min']].astype(int)
#         dict_output_group.update({group: df_group})
#     # final result: {group:df_group}
#     return dict_output_group

# # pack up output of `run_suews`
# # NOT SUITABLE ANYMORE
# def pack_df_output(dict_output):
#     df_raw = pd.DataFrame.from_dict(dict_output).T.applymap(
#         lambda dict: np.concatenate(dict.values())).stack()
#     index = df_raw.index.swaplevel().set_names(['grid', 'tstep'])
#     columns = gen_MultiIndex(dict_output[1][1])
#     values = np.vstack(df_raw.values)
#     df_output = pd.DataFrame(values, index=index, columns=columns)
#     return df_output
#
# def pack_df_state(dict_state):
#     df_raw = pd.DataFrame.from_dict(dict_state).unstack()
#     df_raw = df_raw.map(pd.Series)
#     values = np.vstack(df_raw.values)
#     index = df_raw.index.rename(['tstep', 'grid'])
#     columns = df_raw.iloc[0].index
#     df_state = pd.DataFrame(values, index=index, columns=columns)
#     return df_state
#
# # DEPRECATED: this is slow
# # pack up output of `run_suews`
# def pack_df_output_dep1(dict_output):
#     # # pack all grid and times into index/columns
#     # df_xx = df.from_dict(dict_output, orient='index')
#     # # pack
#     # df_xx1 = df_xx.applymap(lambda s: pd.Series(s)).applymap(df.from_dict)
#     # df_xx2 = pd.concat({grid: pd.concat(
#     #     df_xx1[grid].to_dict()).unstack().dropna(axis=1)
#     #     for grid in df_xx1.columns})
#     # # drop redundant levels
#     # df_xx2.columns = df_xx2.columns.droplevel()
#     # # regroup by `grid`
#     # df_xx2.index.names = ['grid', 'time']
#     # gb_xx2 = df_xx2.groupby(level='grid')
#     # # merge results of each grid
#     # xx3 = gb_xx2.agg(lambda x: tuple(x.values)).applymap(np.array)
#
#     # repack for later concatenation
#     dict_xx4 = ({k: pd.concat(v)
#                  for k, v
#                  in pack_df_grid(dict_output).applymap(df).to_dict().items()})
#
#     # concatenation across model groups
#     res_concat = []
#     for group_x in (x for x in dict_xx4.keys() if not x == 'datetimeline'):
#         # print group_x
#         xx5 = pd.concat((dict_xx4['datetimeline'], dict_xx4[group_x]), axis=1)
#         xx5.columns = gen_group_cols(group_x)
#         res_concat.append(xx5)
#
#     # concatenation across model groups
#     df_output = pd.concat(res_concat, axis=1)
#     # add index information
#     df_output.index.names = ['grid', 'tstep']
#
#     return df_output


# # DEPRECATED: this is slow
# # pack up output of `run_suews`
# def pack_df_output_dep(dict_output):
#     # dict_output is the first value returned by `run_suews`
#     df_res_grid = pd.DataFrame(dict_output).T.stack().swaplevel()
#     dict_grid_time = {grid: pack_dict_output_grid(
#         df_res_grid[grid]) for grid in df_res_grid.index.get_level_values(0)}
#     df_grid_group = pd.DataFrame(dict_grid_time).T
#     return df_grid_group


# DEPRECATED run_suews
# def run_suews(df_forcing, dict_init):
#     # initialise dicts for holding results and model states
#     dict_output = {}
#     # dict_state = {}
#     # dict_state_grid = {grid: dict_state
#     #                    for grid, dict_state
#     # in copy.deepcopy(dict_init).items()}
#     # dict_state is used to save model states for later use
#     t_start = df_forcing.index[0]
#     dict_state = {t_start: copy.deepcopy(dict_init)}
#     # temporal loop
#     for tstep in df_forcing.index:
#         # print 'tstep at', tstep
#         # initialise output of tstep:
#         dict_output.update({tstep: {}})
#         # dict_state is used to save model states for later use
#         dict_state.update({tstep + 1: {}})
#         # load met_forcing if the same across all grids:
#         met_forcing_tstep = df_forcing.loc[tstep].to_dict()
#         # met_forcing_tstep[
#         #     'metforcingdata_grid'] = dict_forcing['metforcingdata_grid']
#         # met_forcing_tstep = df_forcing.iloc[tstep]
#         # xday = met_forcing_tstep['id']
#         # print 'dict_state', dict_state[tstep]
#
#         # spatial loop
#         for grid in dict_state[tstep].keys():
#             dict_state_start = dict_state[tstep][grid]
#             # print 'start', sorted(dict_state_start.keys())
#             # print 'dict_state_start cpanohm', dict_state_start['cpanohm'][0]
#             # print 'dict_state_start xwf', dict_state_start['xwf']
#             # print 'dict_state_start storageheatmethod', dict_state_start['storageheatmethod']
#             # print 'start', dict_state[tstep][grid]['dayofweek'][xday]
#             # mod_config = dict_init[grid]['mod_config']
#             # xx=sp.suews_cal_tstep(
#             #     dict_state_start, met_forcing_tstep, mod_config)
#
#             # calculation at one step:
#             state_end, output_tstep = suews_cal_tstep(
#                 dict_state_start, met_forcing_tstep)
#             # update model state
#             # dict_state_grid[grid].update(state_end)
#             # print 'end', dict_state_grid[grid]['state'][0]
#
#             # update output & model state at tstep for the current grid
#             dict_output[tstep].update({grid: output_tstep})
#             # dict_state[tstep + 1].update({grid: copy.deepcopy(state_end)})
#             # update model state
#             dict_state[tstep + 1].update(
#                 {grid: {k: v for k, v in state_end.items()}})
#
#             # print ''
#             # if tstep==dict_forcing.keys()[-1]:
#             #     print dict_state[tstep][grid]['dayofweek'][43:46]
#
#     return dict_output, dict_state


# # look up properties according to code
# def lookup_code_sub(libName, codeKey, codeValue):
#     # print ''
#     # print 'lookup_code_sub'
#     # print 1,codeName
#     # print 2,codeKey
#     # print 3,codeValue
#     # print ''
#     str_lib = dict_Code2File[libName].replace(
#         '.txt', '').replace('SUEWS', 'lib')
#     str_code = '{:d}'.format(int(np.unique(codeValue)))
#     if codeKey == ':':
#         cmd = '{lib}.loc[{code},:].tolist()'.format(lib=str_lib, code=str_code)
#     else:
#         cmd = '{lib}.loc[{code},{key:{c}^{n}}]'.format(
#             lib=str_lib, code=str_code,
#             key=codeKey, n=len(codeKey) + 2, c='\'')
#     # print cmd
#     res = eval(cmd)
#     return res

# functions for forcing resampling
#
# def EquationOfTime(day):
#     b = (2 * math.pi / 364.0) * (day - 81)
#     return (9.87 * math.sin(2 * b)) - (7.53 * math.cos(b)) - (1.5 * math.sin(b))
#
#
# # r is earth radius vector [astronomical units]
# def GetAberrationCorrection(radius_vector):
#     return -20.4898 / (3600.0 * radius_vector)
#
#
# def GetAltitudeFast(latitude_deg, longitude_deg, utc_datetime):
#     # expect 19 degrees for
#     # solar.GetAltitude(42.364908,-71.112828,datetime.datetime(2007, 2, 18,
#     # 20, 13, 1, 130320))
#
#     day = GetDayOfYear(utc_datetime)
#     declination_rad = math.radians(GetDeclination(day))
#     latitude_rad = math.radians(latitude_deg)
#     hour_angle = GetHourAngle(utc_datetime, longitude_deg)
#
#     first_term = math.cos(
#         latitude_rad) * math.cos(declination_rad) * math.cos(math.radians(hour_angle))
#     second_term = math.sin(latitude_rad) * math.sin(declination_rad)
#     return math.degrees(math.asin(first_term + second_term))
#
#
# def GetCoefficient(jme, constant_array):
#     return sum(
#         [constant_array[i - 1][0] * math.cos(constant_array[i - 1][1] + (constant_array[i - 1][2] * jme)) for i in
#          range(len(constant_array))])
#
#
# def GetDayOfYear(utc_datetime):
#     year_start = datetime(
#         utc_datetime.year, 1, 1, tzinfo=utc_datetime.tzinfo)
#     delta = (utc_datetime - year_start)
#     return delta.days
#
#
# def GetDeclination(day):
#     return 23.45 * math.sin((2 * math.pi / 365.0) * (day - 81))
#
#
# def GetHourAngle(utc_datetime, longitude_deg):
#     solar_time = GetSolarTime(longitude_deg, utc_datetime)
#     return 15 * (12 - solar_time)
#
#
# def GetSolarTime(longitude_deg, utc_datetime):
#     day = GetDayOfYear(utc_datetime)
#     return (((utc_datetime.hour * 60) + utc_datetime.minute + (4 * longitude_deg) + EquationOfTime(day)) / 60)
#
#
# def random_x_N(x, N):
#     list_N = np.zeros(N)
#     pos = random.sample(np.arange(N), x)
#     list_N[pos] = 1
#     return list_N
#
#
# # randomly distribute rainfall in rainAmongN sub-intervals
# def process_rainAmongN(rain, rainAmongN):
#     if rainAmongN <= 3:
#         scale = 3. / rainAmongN
#     else:
#         scale = 1.
#
#     rain_proc = rain.copy()
#     rain_sub = rain_proc[rain_proc > 0]
#     rain_sub_ind = rain_sub.groupby(rain_sub).groups.values()
#     rain_sub_indx = np.array(
#         [x for x in rain_sub_ind if len(x) == 3]).flatten()
#     rain_sub = rain_proc[rain_sub_indx]
#     rain_sub = np.array([scale * random_x_N(rainAmongN, 3) * sub
#                          for sub in rain_sub.values.reshape(-1, 3)])
#     rain_proc[rain_sub_indx] = rain_sub.flatten()
#
#     return rain_proc

# DEPRECATED: Series is slow compared with dict
# high-level wrapper: suews_cal_tstep
# def suews_cal_tstep_df(series_state_start, series_met_forcing_tstep):
#     # use single dict as input for suews_cal_main
#     # series_input_raw = pd.concat(
#     #     [series_state_start.to_dict(),
#     #      series_met_forcing_tstep.rename(series_state_start.name)])
#     dict_input = series_state_start.to_dict()
#     dict_input.update(series_met_forcing_tstep)
#     # print dict_input.keys()
#     # pick only necessary input variables
#     dict_input = {k: dict_input[k] for k in set_var_input}
#
#     # main calculation:
#     res_suews_tstep = sd.suews_cal_main(**dict_input)
#
#     # update state variables
#     # series_state_end = series_input.loc[series_state_start.index]
#     dict_state_end = {k: dict_input[k] for k in series_state_start.index}
#     series_state_end = pd.Series(dict_state_end)
#
#     # pack output
#     # dict_output
#     series_output = pd.Series(data=res_suews_tstep, index=list_var_output)
#     # dict_output = {k: v for k, v in zip(
#     #     list_var_output, res_suews_tstep)}
#
#     return series_state_end, series_output


# def load_SUEWS_MetForcing_df(dir_start):
#     df_forcing = pd.read_table(fileX, delim_whitespace=True,
#                                comment='!',
#                                error_bad_lines=True
#                                # parse_dates={'datetime': [0, 1, 2, 3]},
#                                # keep_date_col=True,
#                                # date_parser=func_parse_date
#                                ).dropna()
#
#     # convert unit
#     df_forcing['press'] = df_forcing['press'] * 10.
#
#     # two datetime's
#     # df_forcing_shift = df_forcing.copy()
#     # df_forcing_shift.loc[:,
#     #                      ['%' + 'iy', 'id', 'it', 'imin']] = (
#     #     df_forcing_shift.loc[
#     #         :, ['%' + 'iy', 'id', 'it', 'imin']
#     #     ].shift(1).fillna(method='backfill'))
#
#     # pack all records of `id` into `all` as required by AnOHM and others
#     # df_grp = df_forcing_shift.groupby('id')
#     df_grp = df_forcing.groupby('id')
#     dict_id_all = {xid: df_grp.get_group(xid)
#                    for xid in df_forcing['id'].unique()}
#     id_all = df_forcing['id'].apply(lambda xid: dict_id_all[xid])
#     df_merged = df_forcing.merge(id_all.to_frame(name='all'),
#                                  left_index=True,
#                                  right_index=True)
#
#     # rename column names to conform with calling function
#     df_merged = df_merged.rename(columns={
#         '%' + 'iy': 'iy',
#         'id': 'id',
#         'it': 'it',
#         'imin': 'imin',
#         'Kdn': 'avkdn',
#         'RH': 'avrh',
#         'Wind': 'avu1',
#         'fcld': 'fcld_obs',
#         'lai_hr': 'lai_obs',
#         'ldown': 'ldown_obs',
#         'rain': 'precip',
#         'press': 'press_hpa',
#         'QH': 'qh_obs',
#         'Q*': 'qn1_obs',
#         'snow': 'snow_obs',
#         'Td': 'temp_c',
#         'all': 'metforcingdata_grid',
#         'xsmd': 'xsmd'})
#
#     # print df_merged.columns
#     # new columns for later use in main calculation
#     df_merged[['iy', 'id', 'it', 'imin']] = df_merged[[
#         'iy', 'id', 'it', 'imin']].astype(np.int64)
#     df_merged['dectime'] = (df_merged['id'] +
#                             (df_merged['it']
#                              + df_merged['imin'] / 60.) / 24.)
#     df_merged['id_prev_t'] = df_merged['id'].shift(1).fillna(method='backfill')
#     df_merged['iy_prev_t'] = df_merged['iy'].shift(1).fillna(method='backfill')
#
#     df_merged['ts5mindata_ir'] = df_merged['temp_c']
#
#     return df_merged

# print dict_var2SiteSelect
# dict_var2SiteSelect = {
#     'lat': 'lat',
#     'lng': 'lng',
#     'timezone': 'Timezone',
#     'alt': 'Alt',
#     'z': 'z',
#     'snowprof':
#         [{'SnowClearingProfWD': ':'}, {'SnowClearingProfWE': ':'}],
#     'ahprof':
#         {'AnthropogenicCode': [
#             {'EnergyUseProfWD': ':'}, {'EnergyUseProfWE': ':'}]},
#     'popprof':
#         {'AnthropogenicCode': [{'PopProfWD': ':'}, {'PopProfWE': ':'}]},
#     'traffprof':
#         {'AnthropogenicCode': [{'TraffProfWD': ':'}, {'TraffProfWE': ':'}]},
#     'humactivity':
#         {'AnthropogenicCode': [
#             {'ActivityProfWD': ':'}, {'ActivityProfWE': ':'}]},
#     'wuprofa':
#         [{'WaterUseProfAutoWD': ':'}, {'WaterUseProfAutoWE': ':'}],
#     'wuprofm':
#         [{'WaterUseProfManuWD': ':'}, {'WaterUseProfManuWE': ':'}],
#     'ah_min': {'AnthropogenicCode': ['AHMin_WD', 'AHMin_WE']},
#     'ah_slope_cooling':
#     {'AnthropogenicCode': ['AHSlope_Cooling_WD', 'AHSlope_Cooling_WE']},
#     'ah_slope_heating':
#     {'AnthropogenicCode': ['AHSlope_Heating_WD', 'AHSlope_Heating_WE']},
#     'alb': [{'Code_Paved': 'AlbedoMax'},
#             {'Code_Bldgs': 'AlbedoMax'},
#             {'Code_EveTr': 'AlbedoMax'},
#             {'Code_DecTr': 'AlbedoMax'},
#             {'Code_Grass': 'AlbedoMax'},
#             {'Code_Bsoil': 'AlbedoMax'},
#             {'Code_Water': 'AlbedoMax'}],
#     'albmax_evetr': {'Code_EveTr': 'AlbedoMax'},
#     'albmax_dectr': {'Code_DecTr': 'AlbedoMax'},
#     'albmax_grass': {'Code_Grass': 'AlbedoMax'},
#     'albmin_evetr': {'Code_EveTr': 'AlbedoMin'},
#     'albmin_dectr': {'Code_DecTr': 'AlbedoMin'},
#     'albmin_grass': {'Code_Grass': 'AlbedoMin'},
#     'alpha_bioco2':
#     [{'Code_EveTr': {'BiogenCO2Code': 'alpha'}},
#             {'Code_DecTr': {'BiogenCO2Code': 'alpha'}},
#             {'Code_Grass': {'BiogenCO2Code': 'alpha'}}],
#     'alpha_enh_bioco2':
#     [{'Code_EveTr': {'BiogenCO2Code': 'alpha_enh'}},
#             {'Code_DecTr': {'BiogenCO2Code': 'alpha_enh'}},
#             {'Code_Grass': {'BiogenCO2Code': 'alpha_enh'}}],
#     'alt': 'Alt',
#     'flowchange': 'FlowChange',
#     'baset': [{'Code_EveTr': 'BaseT'},
#               {'Code_DecTr': 'BaseT'},
#               {'Code_Grass': 'BaseT'}],
#     'basete': [{'Code_EveTr': 'BaseTe'},
#                {'Code_DecTr': 'BaseTe'},
#                {'Code_Grass': 'BaseTe'}],
#     'basethdd': {'AnthropogenicCode': 'BaseTHDD'},
#     'beta_bioco2':
#     [{'Code_EveTr': {'BiogenCO2Code': 'beta'}},
#             {'Code_DecTr': {'BiogenCO2Code': 'beta'}},
#             {'Code_Grass': {'BiogenCO2Code': 'beta'}}],
#     'beta_enh_bioco2':
#     [{'Code_EveTr': {'BiogenCO2Code': 'beta_enh'}},
#             {'Code_DecTr': {'BiogenCO2Code': 'beta_enh'}},
#             {'Code_Grass': {'BiogenCO2Code': 'beta_enh'}}],
#     'bldgh': 'H_Bldgs',
#     'capmax_dec': {'Code_DecTr': 'StorageMax'},
#     'capmin_dec': {'Code_DecTr': 'StorageMin'},
#     'chanohm': [{'Code_Paved': 'AnOHM_Ch'},
#                 {'Code_Bldgs': 'AnOHM_Ch'},
#                 {'Code_EveTr': 'AnOHM_Ch'},
#                 {'Code_DecTr': 'AnOHM_Ch'},
#                 {'Code_Grass': 'AnOHM_Ch'},
#                 {'Code_Bsoil': 'AnOHM_Ch'},
#                 {'Code_Water': 'AnOHM_Ch'}],
#     'cpanohm': [{'Code_Paved': 'AnOHM_Cp'},
#                 {'Code_Bldgs': 'AnOHM_Cp'},
#                 {'Code_EveTr': 'AnOHM_Cp'},
#                 {'Code_DecTr': 'AnOHM_Cp'},
#                 {'Code_Grass': 'AnOHM_Cp'},
#                 {'Code_Bsoil': 'AnOHM_Cp'},
#                 {'Code_Water': 'AnOHM_Cp'}],
#     'crwmax': {'SnowCode': 'CRWMax'},
#     'crwmin': {'SnowCode': 'CRWMin'},
#     'daywat': {'IrrigationCode':
#                ['DayWat(1)',
#                 'DayWat(2)',
#                 'DayWat(3)',
#                 'DayWat(4)',
#                 'DayWat(5)',
#                 'DayWat(6)',
#                 'DayWat(7)']},
#     'daywatper': {'IrrigationCode':
#                   ['DayWatPer(1)',
#                    'DayWatPer(2)',
#                    'DayWatPer(3)',
#                    'DayWatPer(4)',
#                    'DayWatPer(5)',
#                    'DayWatPer(6)',
#                    'DayWatPer(7)']},
#     'dectreeh': 'H_DecTr',
#     'drainrt': 'LUMPS_DrRate',
#     'ef_umolco2perj': {'AnthropogenicCode': 'EF_umolCO2perJ'},
#     'emis': [{'Code_Paved': 'Emissivity'},
#              {'Code_Bldgs': 'Emissivity'},
#              {'Code_EveTr': 'Emissivity'},
#              {'Code_DecTr': 'Emissivity'},
#              {'Code_Grass': 'Emissivity'},
#              {'Code_Bsoil': 'Emissivity'},
#              {'Code_Water': 'Emissivity'}],
#     'enef_v_jkm': {'AnthropogenicCode': 'EnEF_v_Jkm'},
#     'evetreeh': 'H_EveTr',
#     'faibldg': 'FAI_Bldgs',
#     'faidectree': 'FAI_DecTr',
#     'faievetree': 'FAI_EveTr',
#     'faut': {'IrrigationCode': 'Faut'},
#     'fcef_v_kgkm': {'AnthropogenicCode': 'FcEF_v_kgkm'},
#     'frfossilfuel_heat': {'AnthropogenicCode': 'FrFossilFuel_Heat'},
#     'frfossilfuel_nonheat': {'AnthropogenicCode': 'FrFossilFuel_NonHeat'},
#     'g1': {'CondCode': 'G1'},
#     'g2': {'CondCode': 'G2'},
#     'g3': {'CondCode': 'G3'},
#     'g4': {'CondCode': 'G4'},
#     'g5': {'CondCode': 'G5'},
#     'g6': {'CondCode': 'G6'},
#     'gddfull':
#     [{'Code_EveTr': 'GDDFull'},
#             {'Code_DecTr': 'GDDFull'},
#             {'Code_Grass': 'GDDFull'}],
#     'gsmodel': {'CondCode': 'gsModel'},
#     'ie_a': {'IrrigationCode': ['Ie_a1', 'Ie_a2', 'Ie_a3']},
#     'ie_end': {'IrrigationCode': 'Ie_end'},
#     'ie_m': {'IrrigationCode': ['Ie_m1', 'Ie_m2', 'Ie_m3']},
#     'ie_start': {'IrrigationCode': 'Ie_start'},
#     'internalwateruse_h': {'IrrigationCode': 'InternalWaterUse'},
#     'irrfracconif': 'IrrFr_EveTr',
#     'irrfracdecid': 'IrrFr_DecTr',
#     'irrfracgrass': 'IrrFr_Grass',
#     'kkanohm': [{'Code_Paved': 'AnOHM_Kk'},
#                 {'Code_Bldgs': 'AnOHM_Kk'},
#                 {'Code_EveTr': 'AnOHM_Kk'},
#                 {'Code_DecTr': 'AnOHM_Kk'},
#                 {'Code_Grass': 'AnOHM_Kk'},
#                 {'Code_Bsoil': 'AnOHM_Kk'},
#                 {'Code_Water': 'AnOHM_Kk'}],
#     'kmax': {'CondCode': 'Kmax'},
#     'laimax': [{'Code_EveTr': 'LAIMax'},
#                {'Code_DecTr': 'LAIMax'},
#                {'Code_Grass': 'LAIMax'}],
#     'laimin': [{'Code_EveTr': 'LAIMin'},
#                {'Code_DecTr': 'LAIMin'},
#                {'Code_Grass': 'LAIMin'}],
#     #  'lai_obs': '',
#     'laipower': [{'Code_EveTr': ['LeafGrowthPower1', 'LeafGrowthPower2',
#                                  'LeafOffPower1', 'LeafOffPower2']},
#                  {'Code_DecTr': ['LeafGrowthPower1', 'LeafGrowthPower2',
#                                  'LeafOffPower1', 'LeafOffPower2']},
#                  {'Code_Grass': ['LeafGrowthPower1', 'LeafGrowthPower2',
#                                  'LeafOffPower1', 'LeafOffPower2']}],
#     'laitype': [{'Code_EveTr': 'LAIEq'},
#                 {'Code_DecTr': 'LAIEq'},
#                 {'Code_Grass': 'LAIEq'}],
#     'lat': 'lat',
#     'lng': 'lng',
#     'maxconductance': [{'Code_EveTr': 'MaxConductance'},
#                        {'Code_DecTr': 'MaxConductance'},
#                        {'Code_Grass': 'MaxConductance'}],
#     'maxqfmetab': {'AnthropogenicCode': 'MaxQFMetab'},
#     'minqfmetab': {'AnthropogenicCode': 'MinQFMetab'},
#     'min_res_bioco2': [{'Code_EveTr': {'BiogenCO2Code': 'min_respi'}},
#                        {'Code_DecTr': {'BiogenCO2Code': 'min_respi'}},
#                        {'Code_Grass': {'BiogenCO2Code': 'min_respi'}}],
#     'narp_emis_snow': {'SnowCode': 'Emissivity'},
#     'narp_trans_site': 'NARP_Trans',
#     'ohm_coef':
#     [{'Code_Paved': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                      {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                      {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                      {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
#             {'Code_Bldgs': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
#             {'Code_EveTr': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
#             {'Code_DecTr': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
#             {'Code_Grass': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
#             {'Code_Bsoil': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
#             {'Code_Water': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                             {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]},
#             {'SnowCode': [{'OHMCode_SummerWet': ['a1', 'a2', 'a3']},
#                           {'OHMCode_SummerDry': ['a1', 'a2', 'a3']},
#                           {'OHMCode_WinterWet': ['a1', 'a2', 'a3']},
#                           {'OHMCode_WinterDry': ['a1', 'a2', 'a3']}]}],
#     'ohm_threshsw':
#     [{'Code_Paved': 'OHMThresh_SW'},
#             {'Code_Bldgs': 'OHMThresh_SW'},
#             {'Code_EveTr': 'OHMThresh_SW'},
#             {'Code_DecTr': 'OHMThresh_SW'},
#             {'Code_Grass': 'OHMThresh_SW'},
#             {'Code_Bsoil': 'OHMThresh_SW'},
#             {'Code_Water': 'OHMThresh_SW'},
#             {'SnowCode': 'OHMThresh_SW'}],
#     'ohm_threshwd':
#     [{'Code_Paved': 'OHMThresh_WD'},
#             {'Code_Bldgs': 'OHMThresh_WD'},
#             {'Code_EveTr': 'OHMThresh_WD'},
#             {'Code_DecTr': 'OHMThresh_WD'},
#             {'Code_Grass': 'OHMThresh_WD'},
#             {'Code_Bsoil': 'OHMThresh_WD'},
#             {'Code_Water': 'OHMThresh_WD'},
#             {'SnowCode': 'OHMThresh_WD'}],
#     'pipecapacity': 'PipeCapacity',
#     'popdensdaytime': 'PopDensDay',
#     'popdensnighttime': 'PopDensNight',
#     'pormax_dec': {'Code_DecTr': 'PorosityMax'},
#     'pormin_dec': {'Code_DecTr': 'PorosityMin'},
#     'preciplimit': {'SnowCode': 'PrecipLimSnow'},
#     'preciplimitalb': {'SnowCode': 'PrecipLimAlb'},
#     'qf0_beu': ['QF0_BEU_WD', 'QF0_BEU_WE'],
#     'qf_a': {'AnthropogenicCode': ['QF_A_WD', 'QF_A_WE']},
#     'qf_b': {'AnthropogenicCode': ['QF_B_WD', 'QF_B_WE']},
#     'qf_c': {'AnthropogenicCode': ['QF_C_WD', 'QF_C_WE']},
#     'radmeltfact': {'SnowCode': 'RadMeltFactor'},
#     'raincover': 'LUMPS_Cover',
#     'rainmaxres': 'LUMPS_MaxRes',
#     'resp_a': [{'Code_EveTr': {'BiogenCO2Code': 'resp_a'}},
#                {'Code_DecTr': {'BiogenCO2Code': 'resp_a'}},
#                {'Code_Grass': {'BiogenCO2Code': 'resp_a'}}],
#     'resp_b': [{'Code_EveTr': {'BiogenCO2Code': 'resp_b'}},
#                {'Code_DecTr': {'BiogenCO2Code': 'resp_b'}},
#                {'Code_Grass': {'BiogenCO2Code': 'resp_b'}}],
#     'runofftowater': 'RunoffToWater',
#     's1': {'CondCode': 'S1'},
#     's2': {'CondCode': 'S2'},
#     'sathydraulicconduct':
#     [{'Code_Paved': {'SoilTypeCode': 'SatHydraulicCond'}},
#             {'Code_Bldgs': {'SoilTypeCode': 'SatHydraulicCond'}},
#             {'Code_EveTr': {'SoilTypeCode': 'SatHydraulicCond'}},
#             {'Code_DecTr': {'SoilTypeCode': 'SatHydraulicCond'}},
#             {'Code_Grass': {'SoilTypeCode': 'SatHydraulicCond'}},
#             {'Code_Bsoil': {'SoilTypeCode': 'SatHydraulicCond'}},
#             0.],
#     'sddfull': [{'Code_EveTr': 'SDDFull'},
#                 {'Code_DecTr': 'SDDFull'},
#                 {'Code_Grass': 'SDDFull'}],
#     'sfr': ['Fr_Paved',
#             'Fr_Bldgs',
#             'Fr_EveTr',
#             'Fr_DecTr',
#             'Fr_Grass',
#             'Fr_Bsoil',
#             'Fr_Water'],
#     'snowalbmax': {'SnowCode': 'AlbedoMax'},
#     'snowalbmin': {'SnowCode': 'AlbedoMin'},
#     'snowd':
#     [{'Code_Paved': 'SnowLimPatch'},
#             {'Code_Bldgs': 'SnowLimPatch'},
#             {'Code_EveTr': 'SnowLimPatch'},
#             {'Code_DecTr': 'SnowLimPatch'},
#             {'Code_Grass': 'SnowLimPatch'},
#             {'Code_Bsoil': 'SnowLimPatch'},
#             0.],
#     'snowdensmax': {'SnowCode': 'SnowDensMax'},
#     'snowdensmin': {'SnowCode': 'SnowDensMin'},
#     'snowlimbuild': {'Code_Bldgs': 'SnowLimRemove'},
#     'snowlimpaved': {'Code_Paved': 'SnowLimRemove'},
#     'soildepth':
#     [{'Code_Paved': {'SoilTypeCode': 'SoilDepth'}},
#             {'Code_Bldgs': {'SoilTypeCode': 'SoilDepth'}},
#             {'Code_EveTr': {'SoilTypeCode': 'SoilDepth'}},
#             {'Code_DecTr': {'SoilTypeCode': 'SoilDepth'}},
#             {'Code_Grass': {'SoilTypeCode': 'SoilDepth'}},
#             {'Code_Bsoil': {'SoilTypeCode': 'SoilDepth'}},
#             0.],
#     'soilstorecap':
#     [{'Code_Paved': {'SoilTypeCode': 'SoilStoreCap'}},
#             {'Code_Bldgs': {'SoilTypeCode': 'SoilStoreCap'}},
#             {'Code_EveTr': {'SoilTypeCode': 'SoilStoreCap'}},
#             {'Code_DecTr': {'SoilTypeCode': 'SoilStoreCap'}},
#             {'Code_Grass': {'SoilTypeCode': 'SoilStoreCap'}},
#             {'Code_Bsoil': {'SoilTypeCode': 'SoilStoreCap'}},
#             0.],
#     'startdls': 'StartDLS',
#     'enddls': 'EndDLS',
#     'statelimit':
#     [{'Code_Paved': 'StateLimit'},
#             {'Code_Bldgs': 'StateLimit'},
#             {'Code_EveTr': 'StateLimit'},
#             {'Code_DecTr': 'StateLimit'},
#             {'Code_Grass': 'StateLimit'},
#             {'Code_Bsoil': 'StateLimit'},
#             {'Code_Water': 'StateLimit'}],
#     'surf':
#     [[{'Code_Paved': 'StorageMin'},
#       {'Code_Bldgs': 'StorageMin'},
#       {'Code_EveTr': 'StorageMin'},
#       {'Code_DecTr': 'StorageMin'},
#       {'Code_Grass': 'StorageMin'},
#       {'Code_Bsoil': 'StorageMin'},
#       {'Code_Water': 'StorageMin'}],
#             [{'Code_Paved': 'DrainageEq'},
#              {'Code_Bldgs': 'DrainageEq'},
#              {'Code_EveTr': 'DrainageEq'},
#              {'Code_DecTr': 'DrainageEq'},
#              {'Code_Grass': 'DrainageEq'},
#              {'Code_Bsoil': 'DrainageEq'},
#              {'Code_Water': 'DrainageEq'}],
#             [{'Code_Paved': 'DrainageCoef1'},
#              {'Code_Bldgs': 'DrainageCoef1'},
#              {'Code_EveTr': 'DrainageCoef1'},
#              {'Code_DecTr': 'DrainageCoef1'},
#              {'Code_Grass': 'DrainageCoef1'},
#              {'Code_Bsoil': 'DrainageCoef1'},
#              {'Code_Water': 'DrainageCoef1'}],
#             [{'Code_Paved': 'DrainageCoef2'},
#              {'Code_Bldgs': 'DrainageCoef2'},
#              {'Code_EveTr': 'DrainageCoef2'},
#              {'Code_DecTr': 'DrainageCoef2'},
#              {'Code_Grass': 'DrainageCoef2'},
#              {'Code_Bsoil': 'DrainageCoef2'},
#              {'Code_Water': 'DrainageCoef2'}],
#             [{'Code_Paved': 'StorageMax'},
#              {'Code_Bldgs': 'StorageMax'},
#              {'Code_EveTr': 'StorageMax'},
#              {'Code_DecTr': 'StorageMax'},
#              {'Code_Grass': 'StorageMax'},
#              {'Code_Bsoil': 'StorageMax'},
#              {'Code_Water': 'StorageMax'}],
#             [{'Code_Paved': 'StorageMin'},
#              {'Code_Bldgs': 'StorageMin'},
#              {'Code_EveTr': 'StorageMin'},
#              {'Code_DecTr': 'StorageMin'},
#              {'Code_Grass': 'StorageMin'},
#              {'Code_Bsoil': 'StorageMin'},
#              {'Code_Water': 'StorageMin'}]],
#     'surfacearea': 'SurfaceArea',
#     'tau_a': {'SnowCode': 'tau_a'},
#     'tau_f': {'SnowCode': 'tau_f'},
#     'tau_r': {'SnowCode': 'tau_r'},
#     't_critic_cooling':
#     {'AnthropogenicCode': ['TCritic_Cooling_WD', 'TCritic_Cooling_WE']},
#     't_critic_heating':
#         {'AnthropogenicCode': ['TCritic_Heating_WD', 'TCritic_Heating_WE']},
#     'tempmeltfact': {'SnowCode': 'TempMeltFactor'},
#     'th': {'CondCode': 'TH'},
#     'theta_bioco2': [{'Code_EveTr': {'BiogenCO2Code': 'theta'}},
#                      {'Code_DecTr': {'BiogenCO2Code': 'theta'}},
#                      {'Code_Grass': {'BiogenCO2Code': 'theta'}}],
#     'timezone': 'Timezone',
#     'tl': {'CondCode': 'TL'},
#     'trafficrate': ['TrafficRate_WD', 'TrafficRate_WE'],
#     'trafficunits': {'AnthropogenicCode': 'TrafficUnits'},
#     'waterdepth': {'Code_Water': 'WaterDepth'},
#     'waterdist':
#     [[{'WithinGridPavedCode': 'ToPaved'},
#       {'WithinGridBldgsCode': 'ToPaved'},
#       {'WithinGridEveTrCode': 'ToPaved'},
#       {'WithinGridDecTrCode': 'ToPaved'},
#       {'WithinGridGrassCode': 'ToPaved'},
#       {'WithinGridUnmanBSoilCode': 'ToPaved'}],
#             [{'WithinGridPavedCode': 'ToBldgs'},
#              {'WithinGridBldgsCode': 'ToBldgs'},
#              {'WithinGridEveTrCode': 'ToBldgs'},
#              {'WithinGridDecTrCode': 'ToBldgs'},
#              {'WithinGridGrassCode': 'ToBldgs'},
#              {'WithinGridUnmanBSoilCode': 'ToBldgs'}],
#             [{'WithinGridPavedCode': 'ToEveTr'},
#              {'WithinGridBldgsCode': 'ToEveTr'},
#              {'WithinGridEveTrCode': 'ToEveTr'},
#              {'WithinGridDecTrCode': 'ToEveTr'},
#              {'WithinGridGrassCode': 'ToEveTr'},
#              {'WithinGridUnmanBSoilCode': 'ToEveTr'}],
#             [{'WithinGridPavedCode': 'ToDecTr'},
#              {'WithinGridBldgsCode': 'ToDecTr'},
#              {'WithinGridEveTrCode': 'ToDecTr'},
#              {'WithinGridDecTrCode': 'ToDecTr'},
#              {'WithinGridGrassCode': 'ToDecTr'},
#              {'WithinGridUnmanBSoilCode': 'ToDecTr'}],
#             [{'WithinGridPavedCode': 'ToGrass'},
#              {'WithinGridBldgsCode': 'ToGrass'},
#              {'WithinGridEveTrCode': 'ToGrass'},
#              {'WithinGridDecTrCode': 'ToGrass'},
#              {'WithinGridGrassCode': 'ToGrass'},
#              {'WithinGridUnmanBSoilCode': 'ToGrass'}],
#             [{'WithinGridPavedCode': 'ToBSoil'},
#              {'WithinGridBldgsCode': 'ToBSoil'},
#              {'WithinGridEveTrCode': 'ToBSoil'},
#              {'WithinGridDecTrCode': 'ToBSoil'},
#              {'WithinGridGrassCode': 'ToBSoil'},
#              {'WithinGridUnmanBSoilCode': 'ToBSoil'}],
#             [{'WithinGridPavedCode': 'ToWater'},
#              {'WithinGridBldgsCode': 'ToWater'},
#              {'WithinGridEveTrCode': 'ToWater'},
#              {'WithinGridDecTrCode': 'ToWater'},
#              {'WithinGridGrassCode': 'ToWater'},
#              {'WithinGridUnmanBSoilCode': 'ToWater'}],
#             # the last surface type is tricky: needs to determine which goes in
#             # if ToRunoff !=0, use ToRunoff, otherwise use ToSoilStore
#             [{'WithinGridPavedCode': ['ToRunoff', 'ToSoilStore']},
#              {'WithinGridBldgsCode': ['ToRunoff', 'ToSoilStore']},
#              {'WithinGridEveTrCode': ['ToRunoff', 'ToSoilStore']},
#              {'WithinGridDecTrCode': ['ToRunoff', 'ToSoilStore']},
#              {'WithinGridGrassCode': ['ToRunoff', 'ToSoilStore']},
#              {'WithinGridUnmanBSoilCode': ['ToRunoff', 'ToSoilStore']}]
#      ],
#     'wetthresh':
#     [{'Code_Paved': 'WetThreshold'},
#             {'Code_Bldgs': 'WetThreshold'},
#             {'Code_EveTr': 'WetThreshold'},
#             {'Code_DecTr': 'WetThreshold'},
#             {'Code_Grass': 'WetThreshold'},
#             {'Code_Bsoil': 'WetThreshold'},
#             {'Code_Water': 'WetThreshold'}],
#     'year': 'Year',
#     'z': 'z',
#     'z0': 'z0',
#     'zd': 'zd'}


# def proc_met_forcing(df_met_forcing, step_count):
#     met_forcing_tstep = df_met_forcing.iloc[step_count].to_dict()
#     id_x = met_forcing_tstep['id']
#     df_grp = df_met_forcing.groupby('id')
#     all_id = df_grp.get_group(id_x)
#     met_forcing_tstep.update({'all': np.array(all_id.values, order='F')})
#     return met_forcing_tstep


# # a recursive function to retrieve value based on key sequences
# def lookup_KeySeq(indexKey, subKey, indexCode):
#     # print indexKey, subKey, indexCode
#     if type(subKey) is float:
#         res = subKey
#     elif type(subKey) is unicode:
#         res = lookup_code_sub(indexKey, subKey, indexCode)
#     elif type(subKey) is dict:
#         indexKeyX, subKeyX = subKey.items()[0]
#         indexCodeX = lookup_code_sub(indexKey, indexKeyX, indexCode)
#         res = lookup_KeySeq(indexKeyX, subKeyX, indexCodeX)
#     elif type(subKey) is list:
#         res = []
#         for subKeyX in subKey:
#             indexCodeX = indexCode
#             resX = lookup_KeySeq(indexKey, subKeyX, indexCodeX)
#             res.append(resX)
#     # final result
#     return res

# def load_SUEWS_MetForcing_dict(fileX):
#     rawdata_df = load_SUEWS_MetForcing_df_raw(fileX)
#     # dict_met_forcing = rawdata_df.T.to_dict()
#     dict_met_forcing = rawdata_df.to_dict('index')
#     # dict_met_forcing.update(
#     #     {'metforcingdata_grid': np.array(rawdata_df.values,
#     #                                      dtype=np.float, order='F')})
#     return dict_met_forcing

# # resample input foring to tstep required by model
# def resample_forcing(
#         data_raw, tstep_in, tstep_mod, lat, lon, alt, timezone, kdownzen):
#     # reset index as timestamps
#     data_raw.index = data_raw.loc[:, ['iy', 'id', 'it', 'imin']].apply(
#         func_parse_date_row, 1)
#     # shift by half-tstep_in to generate a time series with instantaneous
#     # values
#     data_raw_shift = data_raw.copy().shift(-tstep_in / 2, freq='S')
#
#     # downscale input data to desired time step
#     data_raw_tstep = data_raw_shift.resample(
#         '{tstep}S'.format(tstep=tstep_mod)).interpolate(
#         method='polynomial', order=1).rolling(
#         window=2, center=False).mean()
#
#     # reindex data_tstep to valid range
#     ix = pd.date_range(
#         data_raw.index[0] - timedelta(seconds=tstep_in - tstep_mod),
#         data_raw.index[-1],
#         freq='{tstep}S'.format(tstep=tstep_mod))
#     data_tstep = data_raw_tstep.copy().reindex(
#         index=ix).bfill().ffill().dropna()
#
#     # adjust solar radiation by zenith correction and total amount distribution
#     if kdownzen == 1:
#         data_tstep["avkdn"] = resample_kdn(
#             data_tstep["avkdn"], tstep_mod, timezone, lat, lon, alt)
#
#     # correct rainfall
#     data_tstep['precip'] = resample_precip(
#         data_raw['precip'], tstep_mod, tstep_in)
#
#     # correct temporal information
#     data_tstep['iy'] = data_tstep.index.year
#     data_tstep['id'] = data_tstep.index.dayofyear
#     data_tstep['it'] = data_tstep.index.hour
#     data_tstep['imin'] = data_tstep.index.minute
#
#     # reset index with numbers
#     data_tstep_out = data_tstep.copy().reset_index(drop=True)
#
#     return data_tstep_out


# def load_SUEWS_MetForcing_df_resample(
#         fileX, tstep_in, tstep_mod, lat, lon, alt, timezone, kdownzen):
#     # load raw data
#     df_forcing_met = pd.read_table(fileX, delim_whitespace=True,
#                                    comment='!',
#                                    error_bad_lines=True
#                                    # parse_dates={'datetime': [0, 1, 2, 3]},
#                                    # keep_date_col=True,
#                                    # date_parser=func_parse_date
#                                    ).dropna()
#
#     # convert unit
#     df_forcing_met['press'] = df_forcing_met['press'] * 10.
#
#     # rename column names to conform with calling function
#     df_forcing_met = df_forcing_met.rename(columns={
#         '%' + 'iy': 'iy',
#         'id': 'id',
#         'it': 'it',
#         'imin': 'imin',
#         'Kdn': 'avkdn',
#         'RH': 'avrh',
#         'Wind': 'avu1',
#         'fcld': 'fcld_obs',
#         'lai_hr': 'lai_obs',
#         'ldown': 'ldown_obs',
#         'rain': 'precip',
#         'press': 'press_hpa',
#         'QH': 'qh_obs',
#         'Q*': 'qn1_obs',
#         'snow': 'snow_obs',
#         'Td': 'temp_c',
#         # 'all': 'metforcingdata_grid',
#         'xsmd': 'xsmd'})
#
#     # resample from tstep_in to tstep
#     # met forcing:
#     df_forcing_met_tstep = resample_forcing_met(
#         df_forcing_met, tstep_in, tstep_mod, lat, lon, alt, timezone, kdownzen)
#     # ESTM surface/inner air temp forcing:
#     # df_forcing_estm_tstep = resample_forcing_estm(
#     #     df_forcing_estm, tstep_in, tstep_mod)
#     # merge forcing (met and ESTM) datasets
#     # df_forcing_tstep = df_forcing_met_tstep.merge(
#     #     df_forcing_estm_tstep,
#     #     left_on=['iy', 'id', 'it', 'imin'],
#     #     right_on=['iy', 'id', 'it', 'imin'])
#
#     # pack all records of `id` into `metforcingdata_grid` for AnOHM and others
#     df_forcing_tstep = df_forcing_met_tstep
#     df_grp = df_forcing_tstep.groupby('id')
#     dict_id_all = {xid: df_grp.get_group(xid)
#                    for xid in df_forcing_tstep['id'].unique()}
#     id_all = df_forcing_tstep['id'].apply(lambda xid: dict_id_all[xid])
#     df_merged = df_forcing_tstep.merge(id_all.to_frame(name='metforcingdata_grid'),
#                                        left_index=True,
#                                        right_index=True)
#
#     # print df_merged.columns
#     # new columns for later use in main calculation
#     df_merged[['iy', 'id', 'it', 'imin']] = df_merged[[
#         'iy', 'id', 'it', 'imin']].astype(np.int64)
#     df_merged['dectime'] = (df_merged['id'] +
#                             (df_merged['it']
#                              + df_merged['imin'] / 60.) / 24.)
#     df_merged['id_prev_t'] = df_merged['id'].shift(1).fillna(method='backfill')
#     df_merged['iy_prev_t'] = df_merged['iy'].shift(1).fillna(method='backfill')
#
#     df_merged['ts5mindata_ir'] = df_merged['temp_c']
#
#     return df_merged
#
#
# def load_SUEWS_MetForcing_df_raw(fileX):
#     df_forcing_met = pd.read_table(fileX, delim_whitespace=True,
#                                    comment='!',
#                                    error_bad_lines=True
#                                    # parse_dates={'datetime': [0, 1, 2, 3]},
#                                    # keep_date_col=True,
#                                    # date_parser=func_parse_date
#                                    ).dropna()
#
#     # convert unit
#     df_forcing_met['press'] = df_forcing_met['press'] * 10.
#
#     # two datetime's
#     # df_forcing_shift = df_forcing_met.copy()
#     # df_forcing_shift.loc[:,
#     #                      ['%' + 'iy', 'id', 'it', 'imin']] = (
#     #     df_forcing_shift.loc[
#     #         :, ['%' + 'iy', 'id', 'it', 'imin']
#     #     ].shift(1).fillna(method='backfill'))
#
#     # pack all records of `id` into `all` as required by AnOHM and others
#     # df_grp = df_forcing_shift.groupby('id')
#     df_grp = df_forcing_met.groupby('id')
#     dict_id_all = {xid: df_grp.get_group(xid)
#                    for xid in df_forcing_met['id'].unique()}
#     id_all = df_forcing_met['id'].apply(lambda xid: dict_id_all[xid])
#     df_merged = df_forcing_met.merge(id_all.to_frame(name='all'),
#                                      left_index=True,
#                                      right_index=True)
#
#     # rename column names to conform with calling function
#     df_merged = df_merged.rename(columns={
#         '%' + 'iy': 'iy',
#         'id': 'id',
#         'it': 'it',
#         'imin': 'imin',
#         'Kdn': 'avkdn',
#         'RH': 'avrh',
#         'Wind': 'avu1',
#         'fcld': 'fcld_obs',
#         'lai_hr': 'lai_obs',
#         'ldown': 'ldown_obs',
#         'rain': 'precip',
#         'press': 'press_hpa',
#         'QH': 'qh_obs',
#         'Q*': 'qn1_obs',
#         'snow': 'snow_obs',
#         'Td': 'temp_c',
#         'all': 'metforcingdata_grid',
#         'xsmd': 'xsmd'})
#
#     # print df_merged.columns
#     # new columns for later use in main calculation
#     df_merged[['iy', 'id', 'it', 'imin']] = df_merged[[
#         'iy', 'id', 'it', 'imin']].astype(np.int64)
#     df_merged['dectime'] = (df_merged['id'] +
#                             (df_merged['it']
#                              + df_merged['imin'] / 60.) / 24.)
#     df_merged['id_prev_t'] = df_merged['id'].shift(1).fillna(method='backfill')
#     df_merged['iy_prev_t'] = df_merged['iy'].shift(1).fillna(method='backfill')
#
#     df_merged['ts5mindata_ir'] = df_merged['temp_c']
#
#     return df_merged
