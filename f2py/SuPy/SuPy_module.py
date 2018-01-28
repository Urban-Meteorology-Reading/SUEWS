###########################################################################
# SUEWS for Python
# Authors:
# Ting Sun, ting.sun@reading.ac.uk
# History:
# 20 Jan 2018: first release
###########################################################################

# load dependency modules
import os
import numpy as np
import pandas as pd
import glob
import f90nml
import datetime
from SUEWS_driver import suews_driver as sd
from scipy import interpolate
import collections
import copy

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

# dictionary:
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
        # transpoe laipower:
        df_gridSurfaceChar.loc[xgrid, 'laipower'] = np.array(
            df_gridSurfaceChar.loc[xgrid, 'laipower']).T
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
    # load RunControl variables
    lib_RunControl = load_SUEWS_RunControl(
        os.path.join(dir_input, 'runcontrol.nml'))
    dict_RunControl = lib_RunControl.loc[:, 'runcontrol'].to_dict()
    # some constant values
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

    dict_InitCond = {}
    df_InitCond = pd.DataFrame(data=dict_InitCond)
    # dict for temporally varying states
    # load surface charasteristics
    df_gridSurfaceChar = load_SUEWS_SurfaceChar(dir_input)
    for grid in df_gridSurfaceChar.index:
        for var in df_gridSurfaceChar.columns:
            cmd = '{varX}=df_gridSurfaceChar.loc[{gridX},{varX:{c}^{n}}]'.\
                format(gridX=grid, varX=var, n=len(var) + 2, c='\'')
            exec(cmd)

        # state_init: temporally-varying states from timestep to timestep
        dict_StateInit = {'gdd': 0 * np.ones((ndays + 1, 5), order='F'),
                          'hdd': 0 * np.ones((ndays + 5, 6), order='F'),
                          'icefrac': 0 * np.ones(7),
                          'lai': 0. * np.ones((ndays + 5, 3), order='F'),
                          # snow TODO
                          'meltwaterstore': 0 * np.ones(7, order='F'),
                          'numcapita': (popdensdaytime +
                                        popdensnighttime) / 2.,
                          'porosity': pormax_dec * np.ones(ndays + 1,
                                                           order='F'),
                          'qn1_av_store': -999. * np.ones(2 * nsh + 1),
                          'qn1_s_av_store': -999. * np.ones(2 * nsh + 1),
                          'qn1_s_store': -999. * np.ones(nsh),
                          'qn1_store': -999. * np.ones(nsh),
                          'snowalb': 0. * np.ones(7),
                          'snowdens': 0. * np.ones(7),
                          'snowfrac': 0. * np.ones(7),
                          'snowpack': 0. * np.ones(7),
                          'wu_day': 0. * np.ones(((ndays + 1), 9),
                                                 order='F'),
                          'soilmoist': 0.5 * np.ones(7),
                          'state': -999. * np.ones(7),
                          'tair24hr': 273.15 * np.ones(24 * nsh),
                          'dayofweek': np.ones((ndays + 1, 3), order='F',
                                               dtype=int),
                          'albdectr': albmax_dectr * np.ones(ndays + 1,
                                                             order='F'),
                          'albevetr': albmax_evetr * np.ones(ndays + 1,
                                                             order='F'),
                          'albgrass': albmax_grass * np.ones(ndays + 1,
                                                             order='F'),
                          'decidcap': surf[DecidSurf][5 - 1] * np.ones(
            ndays + 1,
            order='F')}

        # dict_StateInit = {k: v.tolist() for k, v in dict_StateInit.items()}

        # dict with static properties:
        # 1. model settings: dict_ModConfig
        # 2. surface properties
        dict_SurfaceChar = df_gridSurfaceChar.loc[grid, :]
        # print 'shape laipower' ,(df_gridSurfaceChar.loc[grid,'laipower'].shape)
        # print 'len laipower' ,len(dict_SurfaceChar.loc['laipower'])
        # combine them as mod_config
        dict_ModConfig_grid = dict_SurfaceChar.to_dict()
        dict_ModConfig_grid.update(dict_ModConfig)
        dict_ModConfig_grid.update(gridiv=grid)
        dict_ModConfig_grid.keys()
        # construct a DataFrame with entries as:
        # {grid:[mod_config,state_init]}
        dict_InitCond.update(
            {grid: {'mod_config': dict_ModConfig_grid,
                    'state_init': dict_StateInit}})
    return dict_InitCond


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
                               error_bad_lines=True,
                               parse_dates={'datetime': [0, 1, 2, 3]},
                               keep_date_col=True,
                               date_parser=func_parse_date
                               ).dropna()
    df_grp = df_forcing.iloc[:, 1:].groupby('id')
    id_all = df_forcing['id'].apply(lambda xid: df_grp.get_group(xid))
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
    # df_merged['datetime'] = pd.to_datetime()
    df_merged[['iy', 'id', 'it', 'imin']] = df_merged[[
        'iy', 'id', 'it', 'imin']].astype(np.int64)
    df_merged['dectime'] = df_merged['id'] + \
        df_merged['it'] / 24. + df_merged['imin'] / 60.
    df_merged['id_prev_t'] = df_merged['id']
    df_merged['iy_prev_t'] = df_merged['iy']
    df_merged['ts5mindata_ir'] = df_merged['temp_c']

    # convert unit
    df_merged['press_hpa'] = df_merged['press_hpa'] * 10.

    return df_merged


def load_SUEWS_MetForcing_dict(fileX):
    rawdata_df = load_SUEWS_MetForcing_df(fileX)
    # cols = rawdata_df.columns.tolist()
    # rows = rawdata_df.index.tolist()
    # vals = rawdata_df.values.tolist()
    #
    # # dict to hold final result
    # rawdata_dict = {}
    # irow = 0
    # for row in rows:
    #     row_dict = {}
    #     icol = 0
    #     for col in cols:
    #         row_dict.update({col: vals[irow][icol]})
    #         icol += 1
    #     rawdata_dict.update({row: row_dict})
    #     irow += 1

    return rawdata_df.T.to_dict()


def proc_met_forcing(df_met_forcing, step_count):
    met_forcing_tstep = df_met_forcing.iloc[step_count].to_dict()
    id_x = met_forcing_tstep['id']
    df_grp = df_met_forcing.groupby('id')
    all_id = df_grp.get_group(id_x)
    met_forcing_tstep.update({'all': all_id.values})
    return met_forcing_tstep


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


# high-level wrapper: suews_cal_tstep
def suews_cal_tstep(state_old, met_forcing_tstep, mod_config):
    # use single dict as input for suews_cal_main
    # dict_input = met_forcing_tstep.to_dict()
    # dict_input = copy.deepcopy(met_forcing_tstep)
    dict_input = met_forcing_tstep.copy()
    dict_input.update(state_old)
    dict_input.update(mod_config)

    # assign state variables:
    list_var_state = state_old.keys()

    # TODO: ts5mindata_ir needs to be read in from Ts files
    # dict_met['ts5mindata_ir'] = np.array(
    #     dict_met['temp_c'], dtype=np.float, order='F')

    date_time = dict_input['datetime']

    # dayofweek:
    dict_input['dayofweek'] = np.array(
        dict_input['dayofweek']).astype(np.int32)
    id = dict_input['id']
    dict_input['dayofweek'][id, 0] = date_time.dayofweek + 1
    dict_input['dayofweek'][id, 1] = date_time.month
    # season: 1 for summer; 0 for winter.
    dict_input['dayofweek'][id, 2] = (
        1 if 3 < dict_input['dayofweek'][id, 1] < 10 else 0)

    # specify dls:
    dict_input['dls'] = (1 if dict_input['startDLS'] <
                         id < dict_input['endDLS'] else 0)

    # remove unnecessary keys in dict_input
    # redundant keys:
    list_var2del = ['ahprof', 'cbluse', 'disaggmethod', 'disaggmethodestm',
                    'endDLS', 'filecode', 'fileinputpath', 'fileoutputpath',
                    'humactivity', 'datetime',
                    'Kdiff', 'Kdir', 'kdownzen', 'keeptstepfilesin',
                    'keeptstepfilesout',
                    'multipleestmfiles', 'multipleinitfiles',
                    'multiplemetfiles',
                    'popprof',
                    'QE', 'Qf', 'Qs', 'raindisaggmethod', 'resolutionfilesin',
                    'resolutionfilesinestm', 'resolutionfilesout',
                    'solweiguse',
                    'startDLS', 'suppresswarnings', 'traffprof', 'Wd',
                    'writeoutoption', 'wuh', 'wuprofa', 'wuprofm']

    # delete them:
    for k in list_var2del:
        dict_input.pop(k, None)

    # main calculation:
    datetimeline, dataoutline, dataoutlinesnow, dataoutlineestm = \
        sd.suews_cal_main(**dict_input)

    # update state variables
    dict_state = {var: dict_input[var] for var in state_old.keys()}

    # pack output
    dict_output = {'datetime': datetimeline,
                   'dataoutlinesuews': dataoutline,
                   'dataoutlinesnow': dataoutlinesnow,
                   'dataoutlineestm': dataoutlineestm}

    return dict_state, dict_output


# main calculation
def run_suews(dict_forcing, dict_init):
    # initialise dicts for holding results and model states
    dict_output = {}
    # dict_state = {}
    dict_state_grid = {grid: sub_dict['state_init']
                       for grid, sub_dict in copy.deepcopy(dict_init).items()}
    # dict_state is used to save model states for later use
    dict_state = {0: copy.deepcopy(dict_state_grid)}
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

        # spatial loop
        for grid in dict_state_grid.keys():
            state_old = dict_state_grid[grid]
            # print 'start', dict_state_grid[grid]['state']
            # print 'start', dict_state[tstep][grid]['state'][0]
            mod_config = dict_init[grid]['mod_config']
            # xx=sp.suews_cal_tstep(
            #     state_old, met_forcing_tstep, mod_config)
            # calculation at one step:
            state_new, output_tstep = suews_cal_tstep(
                state_old, met_forcing_tstep, mod_config)
            # update model state
            dict_state_grid[grid].update(state_new)
            # print 'end', dict_state_grid[grid]['state'][0]

            # update output & model state at tstep for the current grid
            dict_output[tstep].update({grid: output_tstep})
            # dict_state[tstep + 1].update({grid: copy.deepcopy(state_new)})
            dict_state[tstep + 1].update(
                {grid: {k: v.copy() for k, v in state_new.items()}})
            # print 'dict_state', dict_state[tstep][grid]['state'][0],dict_state[tstep + 1][grid]['state'][0]
            # print ''

    return dict_output, dict_state


# pack up output of one grid of all tsteps
def pack_dict_output_grid(df_grid):
    # get variable info as a DataFrame
    var_df = get_output_info_df()
    # merge dicts of all tsteps
    dict_group = collections.defaultdict(list)
    for d in df_grid:
        for k, v in d.iteritems():  # d.items() in Python 3+
            dict_group[k].append(v)
    # pick groups except for `datatime`
    group_out = (group for group in dict_group.keys()
                 if not group == 'datetime')
    # initialise dict for holding packed output
    dict_output_group = {}
    # pack up output of all tsteps into output groups
    for group_x in group_out:
        # get correct group name by cleaning and swapping case
        group = group_x.replace('dataoutline', '')
        group = (
            group if group in var_df.index.levels[0].drop('datetime')
            else group.swapcase())
        header_group = np.apply_along_axis(
            list, 0, var_df.loc[['datetime', group]].index.values)[:, 1]
        df_group = pd.DataFrame(
            np.hstack((dict_group['datetime'], dict_group[group_x])),
            columns=header_group)
        dict_output_group.update({group: df_group})
    # final result: {group:df_group}
    return dict_output_group

# pack up output of `run_suews`
def pack_df_output(dict_output):
    # dict_output is the first value returned by `run_suews`
    df_res_grid = pd.DataFrame(dict_output).T.stack().swaplevel()
    dict_grid_time = {grid: pack_dict_output_grid(
        df_res_grid[grid]) for grid in df_res_grid.index.get_level_values(0)}
    df_grid_group = pd.DataFrame(dict_grid_time).T
    return df_grid_group


# # test part
# dir_input = './input'
# df_InitCond_grid = init_SUEWS_df(dir_input)
# state_old = df_InitCond_grid.loc[1, 'state_init']
# df_gridSurfaceChar = load_SUEWS_SurfaceChar(dir_input)
#
# # load meterological forcing data: met_forcing_array
# # filecode and resolutionfilesin gives the name to locate met forcing file
# # filecode
# # resolutionfilesin
# list_file_MetForcing = glob.glob(os.path.join(
#     'Input', '{}*{}*txt'.format(filecode, resolutionfilesin / 60)))
# df_forcing = load_SUEWS_MetForcing_df(list_file_MetForcing[0])
# met_forcing = df_forcing
# state_old = df_InitCond_grid.loc[1, 'state_init']
# met_forcing_tstep = met_forcing.loc[1].to_dict()
# met_forcing_tstep.update({'all': np.array(met_forcing.iloc[:, 1:],
#                                           dtype=np.float,
#                                           order='F')})
# mod_config = df_InitCond_grid.loc[1, 'mod_config']
#
# out_state, out_res = suews_cal_tstep(state_old, met_forcing_tstep, mod_config)
# print out_state.keys()
# print out_res.keys()
