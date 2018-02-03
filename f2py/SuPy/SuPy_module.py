###########################################################################
# SUEWS for Python
# Authors:
# Ting Sun, ting.sun@reading.ac.uk
# History:
# 20 Jan 2018: first alpha release
# 01 Feb 2018: performance improvement
# 03 Feb 2018: improvement in output processing
###########################################################################

# load dependency modules
import os
import numpy as np
import pandas as pd
import f90nml
from pandas import DataFrame as df
from SUEWS_driver import suews_driver as sd
from scipy import interpolate
import collections
import copy


######################################################################
# get_args_suews can get the interface informaiton
# of the f2py-converted Fortran interface
def get_args_suews():
    # split doc lines for processing
    docLines = np.array(
        sd.suews_cal_main.__doc__.splitlines(),
        dtype=str)

    # get the information of input variables for SUEWS_driver
    posInput = np.where(
        np.logical_or(
            docLines == 'Parameters', docLines == 'Returns'))
    varInputLines = docLines[posInput[0][0] + 2:posInput[0][1] - 1]
    varInputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                             for x in varInputLines])
    dict_InputInfo = {xx[0]: xx[1] for xx in varInputInfo}
    dict_InOutInfo = {xx[0]: xx[1] for xx in varInputInfo if 'in/out' in xx[1]}

    # get the information of output variables for SUEWS_driver
    posOutput = np.where(docLines == 'Returns')
    varOutputLines = docLines[posOutput[0][0] + 2:]
    varOutputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                              for x in varOutputLines])
    dict_OutputInfo = {xx[0]: xx[1] for xx in varOutputInfo}

    # pack in/out results:
    dict_inout_sd = {
        # 'input' and 'output' are dict's that store variable information:
        # 1. intent: e.g., input, in/output
        # 2. dimension: e.g., (366,7)
        'input': dict_InputInfo,
        'output': dict_OutputInfo,
        # 'var_input' and 'var_output' are tuples,
        # that keep the order of arguments as in the Fortran subroutine
        'var_input': tuple(varInputInfo[:, 0]),
        'var_inout': tuple(dict_InOutInfo.keys()),
        'var_output': tuple(varOutputInfo[:, 0])}

    return dict_inout_sd


##############################################################################
# input processor
# 1. surface properties will be retrieved and packed together for later use
# 2. met forcing conditions will splitted into time steps and used to derive
# other information

# descriptive list/dicts for variables
# minimal required input files for configuration:
list_file_input = ['SUEWS_AnthropogenicHeat.txt',
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

# dictionaries:
# links between code in SiteSelect to properties in according tables
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

# variable translation as done in Fortran-SUEWS
dict_var2SiteSelect = {
    'lat': 'lat',
    'lng': 'lng',
    'timezone': 'Timezone',
    'alt': 'Alt',
    'z': 'z',
    'snowprof':
        [{'SnowClearingProfWD': ':'}, {'SnowClearingProfWE': ':'}],
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
    'startdls': 'StartDLS',
    'enddls': 'EndDLS',
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
            # the last surface type is tricky: needs to determine which goes in
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


# load model settings
# load configurations: mod_config
# process RunControl.nml
# this function can handle all SUEWS nml files
def load_SUEWS_nml(xfile):
    df = pd.DataFrame(f90nml.read(xfile))
    return df

# lib_RunControl = load_SUEWS_nml(os.path.join(dir_input, 'runcontrol.nml'))


def load_SUEWS_RunControl(xfile):
    lib_RunControl = load_SUEWS_nml(xfile)
    for var in lib_RunControl.index:
        val = lib_RunControl.loc[var, 'runcontrol']
        if type(val) == str:
            cmd = '{var}={val:{c}^{n}}'.format(
                var=var, val=val, n=len(val) + 2, c='\'')
        else:
            cmd = '{var}={val}'.format(var=var, val=val)
        # print cmd
        # put configuration variables into global namespace
        exec(cmd, globals())
    # return DataFrame containing settings
    return lib_RunControl


# load all tables (xgrid.e., txt files)
def load_SUEWS_table(fileX):
    rawdata = pd.read_table(fileX, delim_whitespace=True,
                            comment='!', error_bad_lines=True,
                            skiprows=1, index_col=0).dropna()
    return rawdata


# load all tables into variables staring with 'lib_' and filename
def load_SUEWS_Vars(dir_input):
    for k, v in dict_libVar2File.iteritems():
        v = os.path.join(dir_input, v)
        cmd = '{var}=load_SUEWS_table({val:{c}^{n}})'.format(
            var=k, val=v, n=len(v) + 2, c='\'')
        # print cmd
        # put configuration variables into global namespace
        exec(cmd, globals())
    # return DataFrame containing settings
    return None


# look up properties according to code
def lookup_code_sub(codeName, codeKey, codeValue):
    str_lib = dict_Code2File[codeName].replace(
        '.txt', '').replace('SUEWS', 'lib')
    str_code = '{:d}'.format(int(codeValue))
    if codeKey == ':':
        cmd = '{lib}.loc[{code},:].tolist()'.format(lib=str_lib, code=str_code)
    else:
        cmd = '{lib}.loc[{code},{key:{c}^{n}}]'.format(
            lib=str_lib, code=str_code,
            key=codeKey, n=len(codeKey) + 2, c='\'')
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


# load surface characteristics
def load_SUEWS_SurfaceChar(dir_input):
    # load RunControl variables
    lib_RunControl = load_SUEWS_RunControl(
        os.path.join(dir_input, 'runcontrol.nml'))
    dict_RunControl = lib_RunControl.loc[:, 'runcontrol'].to_dict()
    tstep = dict_RunControl['tstep']
    # load all libraries
    load_SUEWS_Vars(dir_input)
    # construct a dictionary in the form: {grid:{var:value,...}}
    dict_gridSurfaceChar = {
        grid: {k: lookup_KeySeq(k, v, grid)
               for k, v in dict_var2SiteSelect.iteritems()}
        for grid in lib_SiteSelect.index}
    # convert the above dict to DataFrame
    df_gridSurfaceChar = pd.DataFrame.from_dict(dict_gridSurfaceChar).T
    # empty dict to hold updated values
    dict_x_grid = {}
    # modify some variables to be compliant with SUEWS requirement
    for xgrid in df_gridSurfaceChar.index:
        # transpoe snowprof:
        df_gridSurfaceChar.loc[xgrid, 'snowprof'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'snowprof'], order='F').T

        # transpoe laipower:
        df_gridSurfaceChar.loc[xgrid, 'laipower'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'laipower'], order='F').T
        # print df_gridSurfaceChar.loc[xgrid, 'laipower'].shape

        # select non-zero values for waterdist of water surface:
        x = np.array(df_gridSurfaceChar.loc[xgrid, 'waterdist'][-1])
        df_gridSurfaceChar.loc[xgrid, 'waterdist'][-1] = (
            x[np.nonzero(x)])

        # surf order as F:
        df_gridSurfaceChar.loc[xgrid, 'surf'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'surf'], order='F')

        # convert to np.array
        df_gridSurfaceChar.loc[xgrid, 'alb'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'alb'])
        # print type(df_gridSurfaceChar.loc[xgrid, 'alb'])

        # dict holding updated values that can be converted to DataFrame later
        dict_x = df_gridSurfaceChar.loc[xgrid, :].to_dict()
        # print 'len(dict_x)',len(dict_x['laipower'])

        # profiles:
        t_tstep = np.linspace(0, 24, num=3600 / tstep * 24, endpoint=False)
        list_varTstep = ['wuprofm_tstep',
                         'ahprof_tstep',
                         'popprof_tstep',
                         'traffprof_tstep',
                         'humactivity_tstep',
                         'wuprofa_tstep']
        for var in list_varTstep:
            var0 = var.replace('_tstep', '')
            var0 = eval(
                'np.array(df_gridSurfaceChar.loc[xgrid, {var:{c}^{n}}]).T'.
                format(var=var0, n=len(var0) + 2, c='\''))
            var0 = np.vstack((var0, var0))
            # interpolator:
            f = interpolate.interp1d(np.arange(0, 48), var0, axis=0)
            cmd = 'dict_x.update({var}=f(t_tstep).tolist())'.format(
                var=var, n=len(var) + 2, c='\'')
            exec(cmd, locals())
            # print cmd

        # update dict to hold grids
        dict_x_grid.update({xgrid: dict_x})

    # convert to DataFrame
    df_x_grid = pd.DataFrame.from_dict(dict_x_grid).T
    return df_x_grid


# create initial conditions
def init_SUEWS_dict(dir_input):  # return dict
    # initialise dict_state_init
    dict_state_init = {}
    # load RunControl variables
    lib_RunControl = load_SUEWS_RunControl(
        os.path.join(dir_input, 'runcontrol.nml'))
    dict_RunControl = lib_RunControl.loc[:, 'runcontrol'].to_dict()
    # some constant values
    nan = -999.
    ndays = 366
    nsh = 3600 / dict_RunControl['tstep']  # tstep from dict_RunControl
    DecidSurf = 4 - 1

    # mod_config: static properties
    dict_ModConfig = {'aerodynamicresistancemethod': 2,
                      'ity': 2,
                      'laicalcyes': 1,
                      'veg_type': 1,
                      'diagqn': 0,
                      'diagqs': 0}
    dict_ModConfig.update(dict_RunControl)

    # dict for temporally varying states
    # load surface charasteristics
    df_gridSurfaceChar = load_SUEWS_SurfaceChar(dir_input)
    for grid in df_gridSurfaceChar.index:
        # initialise dict_InitCond with default values
        dict_InitCond = {
            'dayssincerain':  int(nan),
            'temp_c0':  nan,
            'leavesoutinitially':  int(nan),
            'gdd_1_0':  nan,
            'gdd_2_0':  nan,
            'laiinitialevetr':  nan,
            'laiinitialdectr':  nan,
            'laiinitialgrass':  nan,
            'albevetr0':  nan,
            'albdectr0':  nan,
            'albgrass0':  nan,
            'decidcap0':  nan,
            'porosity0':  nan,
            'pavedstate':  nan,
            'bldgsstate':  nan,
            'evetrstate':  nan,
            'dectrstate':  nan,
            'grassstate':  nan,
            'bsoilstate':  nan,
            'waterstate':  nan,
            'soilstorepavedstate':  nan,
            'soilstorebldgsstate':  nan,
            'soilstoreevetrstate':  nan,
            'soilstoredectrstate':  nan,
            'soilstoregrassstate':  nan,
            'soilstorebsoilstate':  nan,
            'snowinitially':  int(nan),
            'snowwaterpavedstate':  nan,
            'snowwaterbldgsstate':  nan,
            'snowwaterevetrstate':  nan,
            'snowwaterdectrstate':  nan,
            'snowwatergrassstate':  nan,
            'snowwaterbsoilstate':  nan,
            'snowwaterwaterstate':  nan,
            'snowpackpaved':  nan,
            'snowpackbldgs':  nan,
            'snowpackevetr':  nan,
            'snowpackdectr':  nan,
            'snowpackgrass':  nan,
            'snowpackbsoil':  nan,
            'snowpackwater':  nan,
            'snowfracpaved':  nan,
            'snowfracbldgs':  nan,
            'snowfracevetr':  nan,
            'snowfracdectr':  nan,
            'snowfracgrass':  nan,
            'snowfracbsoil':  nan,
            'snowfracwater':  nan,
            'snowdenspaved':  nan,
            'snowdensbldgs':  nan,
            'snowdensevetr':  nan,
            'snowdensdectr':  nan,
            'snowdensgrass':  nan,
            'snowdensbsoil':  nan,
            'snowdenswater':  nan,
            'snowalb0':  nan
        }
        # load Initial Condition variables from namelist file
        lib_InitCond = load_SUEWS_nml(os.path.join(
            dir_input, 'initialconditions{site}_{year}.nml'.format(
                site=dict_ModConfig['filecode'],
                year=int(df_gridSurfaceChar.loc[grid, 'year']))))
        # update default InitialCond with values set in namelist
        dict_InitCond.update(
            lib_InitCond.loc[:, 'initialconditions'].to_dict())

        # snowflag = 1 if snow-related modules are enabled or 0 otherwise
        snowflag = (
            0 if (dict_RunControl['snowuse'] == 0 or
                  dict_InitCond['snowinitially'] == 0)
            else 1)
        # hdd-related parameters:
        dif_basethdd_temp_c0 = (
            df_gridSurfaceChar.loc[
                grid, 'basethdd'] - dict_InitCond['temp_c0'])
        gamma1 = (1 if dif_basethdd_temp_c0 >= 0 else 0)
        gamma2 = (1 if dif_basethdd_temp_c0 <= 0 else 0)
        hdd1 = gamma1 * dif_basethdd_temp_c0
        hdd2 = gamma2 * (-1 * dif_basethdd_temp_c0)
        # state_init: temporally-varying states from timestep to timestep
        dict_StateInit = {
            # water use patterns: TODO: currently not used
            'wu_day': 0. * np.ones(((ndays + 1), 9),
                                   order='F'),
            'numcapita': df_gridSurfaceChar.loc[
                grid,
                ['popdensdaytime', 'popdensnighttime']].mean(),
            'qn1_av_store': nan * np.ones(2 * nsh + 1),
            'qn1_s_av_store': nan * np.ones(2 * nsh + 1),
            'qn1_s_store': nan * np.ones(nsh),
            'qn1_store': nan * np.ones(nsh),

            # snow:
            'snowalb': dict_InitCond['snowalb0'],
            'snowfallcum': 0.,
            'icefrac': 0.2 * np.ones(7),
            # Initial liquid (melted) water for each surface
            'meltwaterstore': snowflag * np.array(
                [dict_InitCond[var] for var in ['snowwaterpavedstate',
                                                'snowwaterbldgsstate',
                                                'snowwaterevetrstate',
                                                'snowwaterdectrstate',
                                                'snowwatergrassstate',
                                                'snowwaterbsoilstate',
                                                'snowwaterwaterstate']],
                order='F'),
            'snowdens': snowflag * np.array(
                [dict_InitCond[var] for var in ['snowdenspaved',
                                                'snowdensbldgs',
                                                'snowdensevetr',
                                                'snowdensdectr',
                                                'snowdensgrass',
                                                'snowdensbsoil',
                                                'snowdenswater']],
                order='F'),
            'snowfrac': snowflag * np.array(
                [dict_InitCond[var] for var in ['snowfracpaved',
                                                'snowfracbldgs',
                                                'snowfracevetr',
                                                'snowfracdectr',
                                                'snowfracgrass',
                                                'snowfracbsoil',
                                                'snowfracwater']],
                order='F'),
            'snowpack': snowflag * np.array(
                [dict_InitCond[var] for var in ['snowpackpaved',
                                                'snowpackbldgs',
                                                'snowpackevetr',
                                                'snowpackdectr',
                                                'snowpackgrass',
                                                'snowpackbsoil',
                                                'snowpackwater']],
                order='F'),

            # Initial soil stores for each surface (below ground)
            'soilmoist': np.array(
                [dict_InitCond[var] for var in ['soilstorepavedstate',
                                                'soilstorebldgsstate',
                                                'soilstoreevetrstate',
                                                'soilstoredectrstate',
                                                'soilstoregrassstate',
                                                'soilstorebsoilstate']]
                + [0.], order='F'),

            # Initial wetness status of each surface (above ground)
            'state': np.array(
                [dict_InitCond[var] for var in ['pavedstate',
                                                'bldgsstate',
                                                'evetrstate',
                                                'dectrstate',
                                                'grassstate',
                                                'bsoilstate',
                                                'waterstate']],
                order='F'),

            # mean Tair of past 24 hours
            'tair24hr': 273.15 * np.ones(24 * nsh),

            # # day of week information:[day, month, season]
            # 'dayofweek': np.ones((ndays + 1, 3), order='F',
            #                      dtype=int),

            # vegetation related parameters:
            'lai': np.vstack((np.array(
                [dict_InitCond[var] for var in ['laiinitialevetr',
                                                'laiinitialdectr',
                                                'laiinitialgrass']],
                order='F'), np.zeros((ndays + 4, 3),
                                     order='F'))),
            'porosity': df_gridSurfaceChar.loc[
                grid, 'pormax_dec'] * np.ones(ndays + 1,
                                              order='F'),
            'albdectr': df_gridSurfaceChar.loc[
                grid, 'albmax_dectr'] * np.ones(ndays + 1,
                                                order='F'),
            'albevetr': df_gridSurfaceChar.loc[
                grid, 'albmax_evetr'] * np.ones(ndays + 1,
                                                order='F'),
            'albgrass': df_gridSurfaceChar.loc[
                grid, 'albmax_grass'] * np.ones(ndays + 1,
                                                order='F'),
            'decidcap': df_gridSurfaceChar.loc[
                grid, 'surf'][DecidSurf][5 - 1] * np.ones(
                ndays + 1,
                order='F'),

            # growing degree days:
            'gdd': 1. * np.vstack((np.array(
                [dict_InitCond[var] for var in ['gdd_1_0',
                                                'gdd_2_0']] + [90, -90, 0],
                order='F'),
                np.array(ndays * [0, 0, 90, -90, 0],
                         order='F').reshape((ndays, -1)))
            ),

            # heating degree days:
            'hdd': 1. * np.vstack((np.array(
                [3 *
                 [0, 0, dict_InitCond['temp_c0'], 0, 0, 0] + \
                 [hdd1, hdd2, dict_InitCond['temp_c0'], 0, 0, \
                  dict_InitCond['dayssincerain']
                  ]
                 ],
                order='F').reshape((-1, 6)),
                np.zeros((ndays + 1, 6), order='F'))),
        }

        dict_StateInit = {k: np.array(v, order='F')
                          for k, v in dict_StateInit.items()}

        # dict with all properties of one grid:
        # 1. model settings: dict_ModConfig
        # 2. surface properties
        # combine them as dict_grid
        dict_grid = df_gridSurfaceChar.loc[grid, :].to_dict()
        dict_grid.update(dict_ModConfig)
        dict_grid.update(dict_StateInit)
        dict_grid.update(gridiv=grid)
        # filter out unnecessary entries for main calculation
        dict_state_init_grid = {
            k: dict_grid[k] for k in list(
                set(dict_grid.keys()).intersection(
                    set(get_args_suews()['var_input'])))}

        # kepp other entries as model configuration
        dict_mod_cfg = {
            k: dict_grid[k] for k in list(
                set(dict_grid.keys()) - set(get_args_suews()['var_input']))}

        # construct a dict with entries as:
        # {grid: dict_state_init}
        dict_state_init.update({grid: dict_state_init_grid})
    # end grid loop

    return dict_mod_cfg, dict_state_init


def init_SUEWS_df(dir_input):  # return pd.DataFrame
    dict_InitCond = init_SUEWS_dict(dir_input)
    df_InitCond = pd.DataFrame.from_dict(dict_InitCond).T

    return df_InitCond


# create met forcing conditions
def func_parse_date(year, doy, hour, min):
    # dt = datetime.datetime.strptime(
    #     ' '.join([year, doy, hour, min]), '%Y %j %H %M')
    dt = pd.to_datetime(' '.join(
        [str(k) for k in [year, doy, hour, min]]),
        format='%Y %j %H %M')
    return dt


def load_SUEWS_MetForcing_df(fileX):
    df_forcing = pd.read_table(fileX, delim_whitespace=True,
                               comment='!',
                               error_bad_lines=True
                               # parse_dates={'datetime': [0, 1, 2, 3]},
                               # keep_date_col=True,
                               # date_parser=func_parse_date
                               ).dropna()
    df_grp = df_forcing.iloc[:, 1:].groupby('id')
    dict_id_all = {xid: df_grp.get_group(xid)
                   for xid in df_forcing['id'].unique()}
    id_all = df_forcing['id'].apply(lambda xid: dict_id_all[xid])
    df_merged = df_forcing.merge(id_all.to_frame(name='all'),
                                 left_index=True,
                                 right_index=True)

    # rename column names to conform with calling function
    df_merged = df_merged.rename(columns={
        '%' + 'iy': 'iy',
        'id': 'id',
        'it': 'it',
        'imin': 'imin',
        'Kdn': 'avkdn',
        'RH': 'avrh',
        'Wind': 'avu1',
        'fcld': 'fcld_obs',
        'lai_hr': 'lai_obs',
        'ldown': 'ldown_obs',
        'rain': 'precip',
        'press': 'press_hpa',
        'QH': 'qh_obs',
        'Q*': 'qn1_obs',
        'snow': 'snow_obs',
        'Td': 'temp_c',
        'xsmd': 'xsmd',
        'all': 'metforcingdata_grid'})

    # new columns for later use in main calculation
    df_merged[['iy', 'id', 'it', 'imin']] = df_merged[[
        'iy', 'id', 'it', 'imin']].astype(np.int64)
    df_merged['dectime'] = df_merged['id'] + \
        df_merged['it'] / 24. + df_merged['imin'] / 60.
    df_merged['id_prev_t'] = df_merged['id'].shift(1).fillna(method='backfill')
    df_merged['iy_prev_t'] = df_merged['iy'].shift(1).fillna(method='backfill')
    # convert unit
    df_merged['press_hpa'] = df_merged['press_hpa'] * 10.
    # TODO: ts5mindata_ir needs to be read in from Ts files
    df_merged['ts5mindata_ir'] = df_merged['temp_c']

    return df_merged


def load_SUEWS_MetForcing_dict(fileX):
    rawdata_df = load_SUEWS_MetForcing_df(fileX)
    return rawdata_df.T.to_dict()


def proc_met_forcing(df_met_forcing, step_count):
    met_forcing_tstep = df_met_forcing.iloc[step_count].to_dict()
    id_x = met_forcing_tstep['id']
    df_grp = df_met_forcing.groupby('id')
    all_id = df_grp.get_group(id_x)
    met_forcing_tstep.update({'all': all_id.values})
    return met_forcing_tstep

# input processing code end here
##############################################################################


##############################################################################
# main calculation
# 1. calculation code for one time step
# 2. compact wrapper for running a whole simulation


# 1. calculation code for one time step
# store these lists for later use
list_var_input = get_args_suews()['var_input']
list_var_inout = get_args_suews()['var_inout']
list_var_output = get_args_suews()['var_output']


# high-level wrapper: suews_cal_tstep
def suews_cal_tstep(dict_state_start, met_forcing_tstep):
    # use single dict as input for suews_cal_main
    dict_input = met_forcing_tstep.copy()
    dict_input.update(dict_state_start)
    # print 'to del:', set(dict_input.keys()) -
    # set(get_args_suews()['var_input'])
    dict_input = {k: dict_input[k] for k in list_var_input}

    # main calculation:
    res_suews_tstep = sd.suews_cal_main(**dict_input)

    # update state variables
    dict_state_end = {var: (copy.copy(dict_input[var])
                            # only copy those variables changed in the fly
                            # to keep better performance
                            if var in list_var_inout
                            # other varialbes are just linked by reference
                            else dict_input[var])
                      for var in dict_state_start.keys()}

    # pack output
    dict_output = {k: v for k, v in zip(
        list_var_output, res_suews_tstep)}

    return dict_state_end, dict_output


# 2. compact wrapper for running a whole simulation
# main calculation
def run_suews(dict_forcing, dict_init):
    # initialise dicts for holding results and model states
    dict_output = {}
    # dict_state = {}
    # dict_state_grid = {grid: dict_state
    #                    for grid, dict_state
    # in copy.deepcopy(dict_init).items()}
    # dict_state is used to save model states for later use
    dict_state = {0: copy.deepcopy(dict_init)}
    # temporal loop
    for tstep in dict_forcing.keys():
        # print 'tstep at', tstep
        # initialise output of tstep:
        dict_output.update({tstep: {}})
        # dict_state is used to save model states for later use
        dict_state.update({tstep + 1: {}})
        # load met_forcing if the same across all grids:
        met_forcing_tstep = dict_forcing[tstep]
        # met_forcing_tstep = df_forcing.iloc[tstep]
        # xday = met_forcing_tstep['id']
        # print 'dict_state', dict_state[tstep]

        # spatial loop
        for grid in dict_state[tstep].keys():
            dict_state_start = dict_state[tstep][grid]
            # print 'start', sorted(dict_state_start.keys())
            # print 'start', dict_state_start
            # print 'start', dict_state[tstep][grid]['dayofweek'][xday]
            # mod_config = dict_init[grid]['mod_config']
            # xx=sp.suews_cal_tstep(
            #     dict_state_start, met_forcing_tstep, mod_config)

            # calculation at one step:
            state_end, output_tstep = suews_cal_tstep(
                dict_state_start, met_forcing_tstep)
            # update model state
            # dict_state_grid[grid].update(state_end)
            # print 'end', dict_state_grid[grid]['state'][0]

            # update output & model state at tstep for the current grid
            dict_output[tstep].update({grid: output_tstep})
            # dict_state[tstep + 1].update({grid: copy.deepcopy(state_end)})
            # update model state
            dict_state[tstep + 1].update(
                {grid: {k: v for k, v in state_end.items()}})

            # print ''
            # if tstep==dict_forcing.keys()[-1]:
            #     print dict_state[tstep][grid]['dayofweek'][43:46]

    return dict_output, dict_state

# main calculation end here
##############################################################################


##############################################################################
# post-processing part
# get variable information from Fortran
def get_output_info_df():
    size_var_list = sd.output_size()
    var_list_x = [np.array(sd.output_name_n(i))
                  for i in np.arange(size_var_list) + 1]
    var_list = np.apply_along_axis(np.vectorize(np.str.strip), 0, var_list_x)

    var_list = np.array(
        filter(lambda var_info: filter(None, var_info), var_list))
    var_df = pd.DataFrame(var_list, columns=['var', 'group', 'aggm'])
    var_dfm = var_df.set_index(['group', 'var'])
    return var_dfm


# get variable info as a DataFrame
# save `var_df` for later use
var_df = get_output_info_df()

# dict as var_df but keys in lowercase
var_df_lower = {group.lower(): group for group in var_df.index.levels[0]}


# pack up output of one grid of all tsteps
def pack_dict_output_grid(df_grid):
    # merge dicts of all tsteps
    dict_group = collections.defaultdict(list)
    for d in df_grid:
        for k, v in d.iteritems():  # d.items() in Python 3+
            dict_group[k].append(v)
    # pick groups except for `datetimeline`
    group_out = (group for group in dict_group.keys()
                 if not group == 'datetimeline')
    # initialise dict for holding packed output
    dict_output_group = {}
    # group names in lower case
    var_df_lower = {group.lower(): group for group in var_df.index.levels[0]}
    # pack up output of all tsteps into output groups
    for group_x in group_out:
        # get correct group name by cleaning and swapping case
        group = group_x.replace('dataoutline', '').replace('line', '')
        # print group
        group = var_df_lower[group]
        header_group = np.apply_along_axis(
            list, 0, var_df.loc[['datetime', group]].index.values)[:, 1]
        # print 'header_group', header_group
        df_group = pd.DataFrame(
            np.hstack((dict_group['datetimeline'], dict_group[group_x])),
            columns=header_group)
        # df_group[[
        #     'Year', 'DOY', 'Hour', 'Min']] = df_group[[
        #         'Year', 'DOY', 'Hour', 'Min']].astype(int)
        dict_output_group.update({group: df_group})
    # final result: {group:df_group}
    return dict_output_group


# generate index for variables in different model groups
def gen_group_cols(group_x):
    # get correct group name by cleaning and swapping case
    group = group_x.replace('dataoutline', '').replace('line', '')
    # print group
    group = var_df_lower[group]
    header_group = np.apply_along_axis(
        list, 0, var_df.loc[['datetime', group]].index.values)[:, 1]

    # generate MultiIndex if not `datetimeline`
    if not group_x == 'datetimeline':
        index_group = pd.MultiIndex.from_product([[group], header_group],
                                                 names=['group', 'var'],
                                                 sortorder=None)
    else:
        index_group = header_group

    return index_group


# merge_grid: useful for both `dict_output` and `dict_state`
def pack_df_grid(dict_output):
    # pack all grid and times into index/columns
    df_xx = df.from_dict(dict_output, orient='index')
    # pack
    df_xx1 = df_xx.applymap(lambda s: pd.Series(s)).applymap(df.from_dict)
    df_xx2 = pd.concat({grid: pd.concat(
        df_xx1[grid].to_dict()).unstack().dropna(axis=1)
        for grid in df_xx1.columns})
    # drop redundant levels
    df_xx2.columns = df_xx2.columns.droplevel()
    # regroup by `grid`
    df_xx2.index.names = ['grid', 'time']
    gb_xx2 = df_xx2.groupby(level='grid')
    # merge results of each grid
    xx3 = gb_xx2.agg(lambda x: tuple(x.values)).applymap(np.array)

    return xx3


# pack up output of `run_suews`
def pack_df_output(dict_output):
    # # pack all grid and times into index/columns
    # df_xx = df.from_dict(dict_output, orient='index')
    # # pack
    # df_xx1 = df_xx.applymap(lambda s: pd.Series(s)).applymap(df.from_dict)
    # df_xx2 = pd.concat({grid: pd.concat(
    #     df_xx1[grid].to_dict()).unstack().dropna(axis=1)
    #     for grid in df_xx1.columns})
    # # drop redundant levels
    # df_xx2.columns = df_xx2.columns.droplevel()
    # # regroup by `grid`
    # df_xx2.index.names = ['grid', 'time']
    # gb_xx2 = df_xx2.groupby(level='grid')
    # # merge results of each grid
    # xx3 = gb_xx2.agg(lambda x: tuple(x.values)).applymap(np.array)

    # repack for later concatenation
    dict_xx4 = ({k: pd.concat(v)
                 for k, v
                 in pack_df_grid(dict_output).applymap(df).to_dict().items()})

    # concatenation across model groups
    res_concat = []
    for group_x in (x for x in dict_xx4.keys() if not x == 'datetimeline'):
        # print group_x
        xx5 = pd.concat((dict_xx4['datetimeline'], dict_xx4[group_x]), axis=1)
        xx5.columns = gen_group_cols(group_x)
        res_concat.append(xx5)

    # concatenation across model groups
    df_output = pd.concat(res_concat, axis=1)
    # add index information
    df_output.index.names = ['grid', 'tstep']

    return df_output


# DEPRECATED: this is slow
# pack up output of `run_suews`
def pack_df_output_dep(dict_output):
    # TODO: add output levels as in the Fortran version
    # dict_output is the first value returned by `run_suews`
    df_res_grid = pd.DataFrame(dict_output).T.stack().swaplevel()
    dict_grid_time = {grid: pack_dict_output_grid(
        df_res_grid[grid]) for grid in df_res_grid.index.get_level_values(0)}
    df_grid_group = pd.DataFrame(dict_grid_time).T
    return df_grid_group


##############################################################################
# auxiliary functions

# convert pandas structures to Python native structures
def conv2PyData(df_x):
    # convert to dict in a split way to expose data/values
    dict_x = df(df_x).to_dict('split')
    # convert data into native Python `list`
    dict_x['data'] = [np.array(var).tolist()
                      for var in dict_x['data']]
    # convert back to DataFrame and then to dict again, which is Python-native
    dict_x_nat = df(**dict_x).to_dict()

    return dict_x_nat
