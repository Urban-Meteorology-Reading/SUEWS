################################################################################
# WFDEI Interpolator
################################################################################
# Purpose:
# interpolate WFDEI data at 1 h scale so as to prepare input data for SUEWS
################################################################################
# Authors:
# Lingbo Xue, Uni. of Reading, L.Xue@student.reading.ac.uk
# Ting Sun, Uni. of Reading, Ting.Sun@reading.ac.uk
################################################################################
# History:
# 20160610 LX: initial version.
# 20160618 TS: refactored with pandas.
# 20160619 TS & LX: solar related problems fixed.
################################################################################
# To do:
# 1. add ability to interpolate at specified temporal scale.
# 2. adpat this code for WATCH data so as to make this code a generic
# interpolaor.
################################################################################

# preload packages
from netCDF4 import Dataset
import numpy as np
from astral import Astral
import os
import sys
import meteo
from datetime import datetime, timedelta, date
from scipy import interpolate
import time
import csv
import pandas as pd


a = Astral()

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


# generate WFDEI filename
def path_WFDEI(directory, var, year, month):
    if var == "Rainf":
        path = directory + var + '_WFDEI_CRU'
        fn = var + '_WFDEI_CRU_' + str(year) + "%02.d" % month + ".nc"
    else:
        path = directory + var + '_WFDEI'
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
            if var == "Rainf":
                fn = path_WFDEI(directory, var, year, month)
            else:
                fn = path_WFDEI(directory, var, year, month)

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


# generate interpolated hourly dataset for SUEWS
def write_SUEWS_forcing_1h(WFDEI_path, output_path, year_start, year_end, lat, lon):
    print('*************** WFDEI Data Interpolator *************** ')
    print('start year: ' + str(year_start))
    print('end year: ' + str(year_end))
    print('input coordinates: (' + '%.2f' % lat + ', ' + '%.2f' % lon + ')')
    print('output path: ' + os.path.expanduser(output_path))
    print(' ')

    # convert user input coordinates into grid index in WFDEI dataset
    glat, glon = lon_lat_grid(lat, lon)
    xlat, xlon = WFDEI_get_city_index(glat, glon)

    # WFDEI variables
    var_list = ["SWdown", "LWdown", "Rainf", "Tair", "PSurf", "Wind", "Qair"]
    # import WFDEI raw data
    data_raw = {var: read_WFDEI(
        WFDEI_path, var, year_start, year_end, xlat, xlon) for var in var_list}

    # join all variable time series into one DataFrame
    ts_data_raw = [k for k in data_raw.itervalues()]
    ts_data_raw = [xx.set_index('time') for xx in ts_data_raw]
    data_raw_3h = pd.concat(ts_data_raw, axis=1, join='outer')

    # expand over the whole date range at the scale of 1 h
    ix = pd.date_range(data_raw_3h.index[
        0], data_raw_3h.index[-1] + timedelta(hours=3), freq="H")
    data_raw_1h = data_raw_3h.reindex(index=ix).resample('1h').mean()

    # create space for processed data
    data_proc_1h = (data_raw_1h.copy())

    # Take off 30-min so that values represent average over previous hour
    sol_elev = np.array([Astral.solar_elevation(
        a, t - timedelta(minutes=30), lat, lon) for t in data_raw_1h["SWdown"].index])
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
        # refill starting & ending missing values
        data_proc_1h[var][:3] = data_proc_1h[var][2]
        data_proc_1h[var] = data_proc_1h[var].fillna(method='pad')

    # Rainf: evenly distribute over the 3-h period
    data_proc_1h['Rainf'] = data_raw_1h['Rainf'].interpolate(
        method='polynomial', order=0).shift(-2)

    # SWdown:
    # force nocturnal values to zero:
    data_proc_1h["SWdown"] = sol_elev_reset * data_proc_1h["SWdown"]
    # rescale based on 3-hourly values:
    avg_3h_SWdown = data_raw_3h["SWdown"].resample('D').mean()
    avg_1h_SWdown = data_proc_1h["SWdown"].resample('D').mean()
    ratio_SWdown = (avg_3h_SWdown / avg_1h_SWdown).reindex(
        index=ix).resample('1h').mean().fillna(method='pad')
    data_proc_1h["SWdown"] = (
        ratio_SWdown * data_proc_1h['SWdown']).fillna(method='pad')

    # fill the first three hours with value at 3:00 h
    for var in var_List_p3h:
        data_proc_1h[var][:3] = data_proc_1h[var][2]
        data_proc_1h[var][-3:] = data_proc_1h[var][-2]

    # export processed data
    header = ["iy", "id", "it", "imin", "qn", "qh", "qe", "qs", "qf", "U", "RH", "Tair", "pres",
              "rain", "kdown", "snow", "ldown", "fcld", "wuh", "xsmd", "lai", "kdiff", "kdir", "wdir"]
    data_out_1h = pd.DataFrame(index=data_proc_1h.index, columns=header)
    var_out_list = ['SWdown', 'LWdown', 'Rainf', 'Tair', 'PSurf', 'Wind',]
    var_out_zip = np.array(
        [var_out_list, ['kdown', 'ldown', 'rain', 'Tair', 'pres', 'U',]]).T
    # fill in variables:
    for p in var_out_zip:
        data_out_1h[p[1]] = data_proc_1h[p[0]]

    # RH calculation:
    data_out_1h['RH'] = 100 * meteo.mixr2rh(meteo.sh2mixr(
        data_proc_1h['Qair']), data_out_1h['pres'], data_out_1h['Tair'])

    # unit conversion:
    # Tair: K -> degC
    data_out_1h['Tair'] = data_out_1h['Tair'] - 273.15
    # rainfall: kg m-2 -> mm 3 x 60 x 60 s / 1000 kg m-3 * 1000 mm m-1
    data_out_1h['rain'] = data_out_1h['rain'] / 3 * 3 * 60 * 60
    # presure: Pa -> kPa
    data_out_1h['pres'] = data_out_1h['pres'] / 1000

    # process timestamps
    data_out_1h['iy'] = data_proc_1h.index.year
    data_out_1h['id'] = data_proc_1h.index.dayofyear
    data_out_1h['it'] = data_proc_1h.index.hour
    data_out_1h['imin'] = data_proc_1h.index.minute

    # replace nan with -999
    data_out_1h = data_out_1h.fillna(value=-999)[1:]

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

    print('*************** Interpolation succeeded! *************** ')



################################################################################
# running section:
# provide parameters here
################################################################################
# read in data from WFDEI files
input_path='/Users/sunt05/Documents/Data/WFDEI/'
# input_path = '/Volumes/DATA-TS/WFDEI/'
output_path = '~/Downloads'
year_start, year_end = 2012, 2012
lat, lon = 51.51, -0.12
# city = a["London"]
# lat, lon = lon_lat_grid(city.latitude, city.longitude)

start = time.time()
write_SUEWS_forcing_1h(input_path, output_path,
                       year_start, year_end, lat, lon)
end = time.time()
print('time used in processing:' + '%.2f' % (end - start) + ' s')
