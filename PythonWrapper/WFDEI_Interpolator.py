##########################################################################
# WFDEI Interpolator
##########################################################################
# Purpose:
# 1. interpolate WFDEI data at hourly scale;
# 2. prepare input data for SUEWS.
##########################################################################
# Authors:
# Lingbo Xue, Uni. of Reading, L.Xue@student.reading.ac.uk
# Ting Sun, Uni. of Reading, Ting.Sun@reading.ac.uk
##########################################################################
# History:
# 20160610 LX: initial version.
# 20160618 TS: refactored with pandas.
# 20160619 TS & LX: solar related problems fixed.
# 20160706 TS & LX: solar parts replaced with Pysolar.
# 20160707 TS: relative humidity modified to be consistent with SUEWS.
# 20160707 TS: interactive input implemented.
##########################################################################
# To do:
# 1. add ability to interpolate at specified temporal scale.
# 2. adpat this code for WATCH data so as to make this code a generic
# interpolaor.
##########################################################################

# preload packages
from netCDF4 import Dataset
import numpy as np
import os
import sys
from datetime import datetime, timedelta, date
from scipy import interpolate
import time
import pandas as pd
from Pysolar import solar


# determing grid index according to coordinates
def lon_lat_grid(lat, lon):
    lon_deci = lon - int(lon)
    lat_deci = lat - int(lat)

    if lon >= 0:
        if 0 <= lon_deci < 0.5:
            lon = int(lon) + 0.25
        else:
            lon = int(lon) + 0.75
    else:
        if -0.5 < lon_deci <= 0:
            lon = -(-int(lon) + 0.25)
        else:
            lon = -(-int(lon) + 0.75)

    if lat >= 0:
        if 0 <= lat_deci < 0.5:
            lat = int(lat) + 0.25
        else:
            lat = int(lat) + 0.75
    else:
        if -0.5 < lat_deci <= 0:
            lat = -(-int(lat) + 0.25)
        else:
            lat = -(-int(lat) + 0.75)

    return lat, lon


# Get City Index: WATCH
def WATCH_get_city_index(lat, lon):
    nc = Dataset("WFD-land-lat-long-z.nc")
    for i in range(0, 67420):
        if nc.variables['Latitude'][i] == lat and nc.variables['Longitude'][i] == lon:
            index = i
            break
    return index


# Get City Index: WFDEI
def WFDEI_get_city_index(lat, lon):
    with open('WFDEI-land-long-lat-height.txt') as f:
        ls = [line.split() for line in f]
    for i in range(7, len(ls)):
        if float(ls[i][0]) == lon and float(ls[i][1]) == lat:
            return int(ls[i][4]), int(ls[i][3])
            break


# calculate saturation vapour pressure [Pa]
def evpsat(T_degC, p_hPa):
    ''' get saturation pressure [Pa] for a given air temperature [degC] and pressure [hPa]'''
    from numpy import log10

    T_degC = max(T_degC, 0.001, key=lambda x: abs(x))
    p_kPa = p_hPa / 10

    if T_degC > 0:
        # For T > 0 degC (f corrects for the fact that we are not dealing with
        # pure water)
        e_mb_pos = 6.1121 * np.exp(((18.678 - T_degC / 234.5)
                                    * T_degC) / (T_degC + 257.14))
        f_pos = 1.00072 + p_kPa * (3.2e-6 + 5.9e-10 * T_degC ** 2)
        esat_hPa_pos = e_mb_pos * f_pos
        esat_hPa = esat_hPa_pos
    else:
        # For T <= 0 degC
        e_mb_neg = 6.1115 * np.exp(((23.036 + T_degC / 333.7)
                                    * T_degC) / (T_degC + 279.82))
        f_neg = 1.00022 + p_kPa * (3.83e-6 + 6.4e-10 * T_degC ** 2)
        esat_hPa_neg = e_mb_neg * f_neg
        esat_hPa = esat_hPa_neg

    esat_Pa = esat_hPa * 100
    return esat_Pa


# calculate actual vapour pressure [Pa]
def eact(qv, p_kPa):
    # Specific gas constant for dry air  (Rv = R/RMM(dry air)) / J K^-1 kg^-1
    Rd = 287.04
    # Specific gas contstant for water vapour (Rv = R/RMM(H20)) / J K^-1 kg^-1
    Rv = 461.5
    # calculate actual vapour pressure [Pa]
    evp_Pa = Rv * qv * p_kPa * 1000 / (Rd + qv * (Rv - Rd))
    return evp_Pa


# convert specific humidity [kg/kg] to relative humidity [%]
def q2rh(qv, p_kPa, T_degC):
    eact_Pa = eact(qv, p_kPa)
    esat_Pa = evpsat(T_degC, p_kPa * 10)
    rh_pct = 100 * eact_Pa / esat_Pa

    return rh_pct


# vectorize q2rh
def vq2rh(qv, p_kPa, T_degC):
    eact_Pa = eact(qv, p_kPa)
    vevpsat = np.vectorize(evpsat)
    esat_Pa = vevpsat(T_degC, p_kPa * 10)
    rh_pct = 100 * eact_Pa / esat_Pa

    return rh_pct

# functions for calculating RH
# from package of meteo
Mw = 18.0160  # molecular weight of water
Md = 28.9660  # molecular weight of dry air


def esat(T):
    ''' get sateration pressure (units [Pa]) for a given air temperature (units [K])'''
    from numpy import log10
    TK = 273.15
    e1 = 101325.0
    logTTK = log10(T / TK)
    x1 = 10.79586 * (1 - TK / T)
    x2 = 5.02808 * logTTK
    x3 = 1.50474 * 1e-4 * (1. - 10**(-8.29692 * (T / TK - 1)))
    x4 = 0.42873 * 1e-3 * (10**(4.76955 * (1 - TK / T)) - 1) - 2.2195983
    xx = x1 - x2 + x3 + x4
    esat = e1 * 10 ** xx
    return esat


def sh2mixr(qv):
    '''conversion from specific humidity (units [kg/kg]) to mixing ratio (units also [kg/kg])'''
    return qv / (1. - qv)


def mixr2rh(mixr, p, T):
    '''purpose: conversion mixing ratio to relative humidity [kg/kg] (not tested)'''
    return mixr * p / ((mixr + Mw / Md) * esat(T))


def sh2rh(qv, p, T):
    '''conversion from specific humidity (units [kg/kg]) to relative humidity in percentage'''
    return mixr2rh(sh2mixr(qv), p, T)


# generate WFDEI filename
def path_WFDEI(directory, var, year, month):
    if var == "Rainf":
        path = os.path.join(directory, var + '_WFDEI_CRU')
        fn = var + '_WFDEI_CRU_' + str(year) + "%02.d" % month + ".nc"
    else:
        path = os.path.join(directory, var + '_WFDEI')
        fn = var + '_WFDEI_' + str(year) + "%02.d" % month + ".nc"

    path = os.path.join(path, fn)
    return path


# import WFDEI data
def read_WFDEI(directory, var, year_start, year_end, xlat, xlon):
    print('reading in ' + var + ':')
    rawdata = []
    for year in range(year_start, year_end + 1):
        print('     working on ' + str(year) + '...')
        for month in range(1, 13):
            # determine file name to read in
            fn = path_WFDEI(directory, var, year, month)
            print(fn)
            # get WFDEI dataset:
            nc = Dataset(fn)

            # read in data:
            for i in range(0, len(nc.variables[var][:, xlat, xlon])):
                # determine date time string
                date = str(year) + "%02.d" % month + \
                    "%02.d" % (i / 8 + 1) + "%02.d" % (i % 8 * 3)

                # note the staring index in WFDEI data is 1 whereas in python is 0
                # so xlat-1 and xlon-1 are needed.
                rawdata.append(
                    (date, nc.variables[var][:, xlat - 1, xlon - 1][i]))
    # convert to time series
    ts_data = pd.DataFrame(rawdata, columns=['time', var])
    ts_data["time"] = pd.to_datetime(ts_data['time'], format="%Y%m%d%H")
    print(var + ' successfully imported!')
    print('\n')
    return ts_data

# load and rearrange WFDEI data for interpolation


def load_WFDEI_3h(WFDEI_path, year_start, year_end, lat, lon):
    print('*************** WFDEI Data Loader *************** ')
    print('start year: ' + str(year_start))
    print('end year: ' + str(year_end))
    print('input coordinates: (' + '%.2f' % lat + ', ' + '%.2f' % lon + ')')
    # convert user input coordinates into grid index in WFDEI dataset
    glat, glon = lon_lat_grid(lat, lon)
    xlat, xlon = WFDEI_get_city_index(glat, glon)

    # WFDEI variables
    var_list = ["SWdown", "LWdown", "Rainf", "Tair", "PSurf", "Wind", "Qair"]
    # import WFDEI raw data
    data_raw = {var: read_WFDEI(
        WFDEI_path, var, year_start, year_end, xlat, xlon) for var in var_list}
    # join all variable time series into one DataFrame
    try:
        ts_data_raw = [k for k in data_raw.itervalues()]
    except AttributeError as e:
        ts_data_raw = [k for k in data_raw.values()]

    ts_data_raw = [xx.set_index('time') for xx in ts_data_raw]
    data_raw_3h = pd.concat(ts_data_raw, axis=1, join='outer')

    print('*************** WFDEI Data Successfully Loaded *************** ')
    return data_raw_3h


# interpolate 3-hourly raw data to hourly results for SUEWS
def process_SUEWS_forcing_1h(data_raw_3h, lat, lon):
    #     print('*************** WFDEI Data Processor *************** ')
    # expand over the whole date range at the scale of 1 h
    ix = pd.date_range(data_raw_3h.index[
        0], data_raw_3h.index[-1] + timedelta(hours=3), freq="H")
    data_raw_1h = data_raw_3h.reindex(index=ix).resample('1h').mean()

    # create space for processed data
    data_proc_1h = data_raw_1h.copy()

    # Take off 30-min so that values represent average over previous hour
    # sol_elev = np.array([Astral.solar_elevation(
    # a, t - timedelta(minutes=30), lat, lon) for t in
    # data_raw_1h["SWdown"].index])

    sol_elev = np.array([solar.GetAltitudeFast(
        lat, lon, t - timedelta(minutes=30)) for t in data_raw_1h["SWdown"].index])

    sol_elev_reset = np.sin(np.radians(sol_elev.copy()))
    sol_elev_reset[sol_elev > 0] = 1
    sol_elev_reset[sol_elev <= 0] = 0

    # normal interpolation for instantaneous variables:
    # i.e., Tair, Wind, Psurf, Qair.
    # these variables have been processed
    var_List_inst = ["Tair", "PSurf", "Wind", "Qair"]
    for var in var_List_inst:
        data_proc_1h[var] = data_raw_1h[var].interpolate(method='polynomial', order=1).rolling(
            window=2, center=False).mean().fillna(method='pad')
        # fill the first three hours with value at 3:00 h
        data_proc_1h[var][0:4] = data_proc_1h[var][3]

    # normal interpolation for variables averaged over previous 3 hours:
    # i.e., SWdown, LWdown, Rainf
    var_List_p3h = ["SWdown", "LWdown", "Rainf"]

    # SWdown, LWdown:
    for var in var_List_p3h[:-1]:
        # convert to 30-min instantaneous values
        x0 = data_raw_1h[var].resample('30min').mean(
        ).interpolate(method='polynomial', order=1)
        # shift to get delta/6, so that values at
        # [t-1,t,t+1]=xt+delta*[1,3,5]/6
        data_proc_1h[var] = x0.shift(-3)[::2]

    # Rainf: evenly distribute over the 3-h period
    data_proc_1h['Rainf'] = data_raw_1h['Rainf'].interpolate(
        method='polynomial', order=0).shift(-2).fillna(method='pad', limit=2).fillna(value=0)

    # SWdown correction:
    # force nocturnal values to zero:
    data_proc_1h["SWdown"] = sol_elev_reset * data_proc_1h["SWdown"]
    # rescale based on 3-hourly values:
    avg_3h_SWdown = data_raw_3h["SWdown"].resample('D').mean()
    avg_1h_SWdown = data_proc_1h["SWdown"].resample('D').mean()
    ratio_SWdown = (avg_3h_SWdown / avg_1h_SWdown).reindex(
        index=ix).resample('1h').mean().fillna(method='pad')
    data_proc_1h["SWdown"] = (
        ratio_SWdown * data_proc_1h['SWdown']).fillna(method='pad')

    # export processed data
    header = ["iy", "id", "it", "imin", "qn", "qh", "qe", "qs", "qf", "U", "RH", "Tair", "pres",
              "rain", "kdown", "snow", "ldown", "fcld", "wuh", "xsmd", "lai", "kdiff", "kdir", "wdir"]
    data_out_1h = pd.DataFrame(index=data_proc_1h.index, columns=header)
    var_out_list = ['SWdown', 'LWdown', 'Rainf', 'Tair', 'PSurf', 'Wind']
    var_out_zip = np.array(
        [var_out_list, ['kdown', 'ldown', 'rain', 'Tair', 'pres', 'U']]).T

    # refill starting & ending missing values
    var_fill_list = ['SWdown', 'LWdown', 'Tair', 'PSurf', 'Wind', 'Qair']
    for var in var_fill_list:
        data_proc_1h[var][:4] = data_proc_1h[var][4]
        data_proc_1h[var] = data_proc_1h[var].fillna(method='pad')

    # fill in variables:
    for p in var_out_zip:
        data_out_1h[p[1]] = data_proc_1h[p[0]]

    # RH calculation:
    data_out_1h['RH'] = vq2rh(data_proc_1h['Qair'],
                              data_proc_1h['PSurf'] / 1000,
                              data_proc_1h['Tair'] - 273.15)

    # unit conversion:
    # Tair: K -> degC
    data_out_1h['Tair'] -= 273.15
    # rainfall: kg m-2 -> mm  60 x 60 s / 1000 kg m-3 * 1000 mm m-1
    data_out_1h['rain'] *= 60 * 60
    # presure: Pa -> kPa
    data_out_1h['pres'] /= 1000

    # process timestamps
    data_out_1h['iy'] = data_proc_1h.index.year
    data_out_1h['id'] = data_proc_1h.index.dayofyear
    data_out_1h['it'] = data_proc_1h.index.hour
    data_out_1h['imin'] = data_proc_1h.index.minute

    # replace nan with -999
    data_out_1h = data_out_1h.fillna(value=-999)

    print('*************** WFDEI Data Successfully Processed *************** ')

    return data_out_1h


def write_SUEWS_forcing_1h(WFDEI_path, output_path, year_start, year_end, lat, lon):
    # load raw 3-hourly data
    data_raw_3h = load_WFDEI_3h(input_path, year_start, year_end, lat, lon)

    # process raw data to hourly forcings for SUEWS
    data_out_1h = process_SUEWS_forcing_1h(data_raw_3h, lat, lon)

    # output files of each year
    print('output files:')
    for year in range(year_start, year_end + 1):
        data_out_1h_year = data_out_1h[lambda df: (
            df.index - timedelta(minutes=60)).year == year]
        file_output_year = os.path.expanduser(
            os.path.join(output_path, 'WFDEI_' + str(year) + '.txt'))
        data_out_1h_year.to_csv(file_output_year, sep=" ",
                                index=False, float_format='%.4f')
        print(file_output_year)

    print('********* WFDEI Data Processing Successfully Finished *********')


##########################################################################
# running section:
# provide parameters here
##########################################################################
# read in data from WFDEI files

# input_path = '/Users/sunt05/Documents/Data/WFDEI/'
# input_path = '/Volumes/DATA-TS/WFDEI/'

# output_path = '~/Downloads'

# year_start, year_end = 2012, 2012

# lat, lon = 51.51, -0.12  # London
# lat, lon = 51.58, -1.8  # Swindon


# start = time.time()
# write_SUEWS_forcing_1h(input_path, output_path,
#                        year_start, year_end, lat, lon)
# end = time.time()
print('time used in processing:' + '%.2f' % (end - start) + ' s')


# interatively get input parameters
while True:
    # WFDEI path:
    while True:
        input_path = raw_input("Please input the path for WFDEI data: ")
        input_path = os.path.realpath(os.path.expanduser(input_path))
        print input_path
        if os.path.lexists(input_path):
            break
        else:
            print "No such directory. Try again..."

    # output path:
    while True:
        output_path = raw_input("Please input the path for output: ")
        output_path = os.path.realpath(os.path.expanduser(output_path))
        print output_path
        if os.path.lexists(output_path):
            break
        else:
            print "No such directory. Try again..."

    # year range:
    while True:
        year_start = int(raw_input(
            "Please input the start year (YYYY): "))
        year_end = int(raw_input(
            "Please input the end year (YYYY): "))
        print(1979 <= year_start <= year_end <= 2014)
        if 1979 <= year_start <= year_end <= 2014:
            break
        else:
            print "Please input valid years. Try again..."

    # coordinates:
    while True:
        lat = float(raw_input(
            "Please input the latitude (in deg): "))
        lon = float(raw_input(
            "Please input the longitude (in deg): "))
        print(-90 < lat < 90 and -180 < lon < 180)
        if -90 < lat < 90 and -180 < lon < 180:
            break
        else:
            print "Please input valid coordinates. Try again..."

    start = time.time()
    write_SUEWS_forcing_1h(input_path, output_path,
                           year_start, year_end, lat, lon)
    end = time.time()
    print('time used in processing:' + '%.2f' % (end - start) + ' s')

    t = raw_input('Do you want to quit? Y/N')
    if t == 'Y' or t == 'y':
        ftp.quit()
        break
