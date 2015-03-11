__author__ = 'Fredrik Lindberg'

# This class will be used to prepare input met data into UMEP

import numpy as np


def leap_year(yy):
    if (yy % 4) == 0:
        if (yy % 100) == 0:
            if (yy % 400) == 0:
                leapyear = 1
            else:
                leapyear = 0
        else:
            leapyear = 1
    else:
        leapyear = 0

    return leapyear


class SuewsDataProcessing:
    def __init__(self):
        pass

    def tofivemin_v1(self, met_old):

        # Time columns
        iy = met_old[:, 0]
        id = met_old[:, 1]
        it = met_old[:, 2]
        imin = met_old[:, 3]

        # first figure out the time res of input file
        dectime0 = id[0] + it[0] / 24 + imin[0] / (60 * 24)
        dectime1 = id[1] + it[1] / 24 + imin[1] / (60 * 24)
        timeres_old = np.round((dectime1 - dectime0) * (60 * 24))
        howmanyfives = int(timeres_old / 5)

        if howmanyfives < 1:  # put code here to make longer timestep
            notused = []
        else:  # interpolate to five minute
            met_new = np.zeros((met_old.shape[0] * howmanyfives, 24))
            for i in range(0, met_old.shape[0]):
                timestep_before = met_old[i - 1, :]
                timestep_now = met_old[i, :]

                if met_old.shape[0] - 1 == i:
                    timestep_next = timestep_now  #This should be fixed if multiple years are run
                else:
                    timestep_next = met_old[i + 1, :]

                changeinold_before = (timestep_before - timestep_now) / howmanyfives
                changeinold_next = (timestep_next - timestep_now) / howmanyfives

                if timestep_next[13] == -999:
                    rainchange = 0
                else:
                    rainchange = timestep_next[13] / howmanyfives
                if timestep_next[18] == -999:
                    waterchange = 0
                else:
                    waterchange = timestep_next[18] / howmanyfives

                met_temp = np.zeros((howmanyfives, 24))
                min = 0
                for j in range(0, howmanyfives):
                    min += 5
                    #index = i * j + j
                    if j * 5 < 30:
                        if i == 0:
                            met_temp[j, :] = met_old[i, :]
                        else:
                            met_temp[j, :] = timestep_now + changeinold_before * ((howmanyfives / 2) - j)
                    else:
                        if i == met_old.shape[0] - 2:
                            met_temp[j, :] = met_old[i, :]
                        else:
                            met_temp[j, :] = timestep_now + changeinold_next * (j - (howmanyfives / 2))
                    met_temp[j, 13] = rainchange  # rain
                    met_temp[j, 18] = waterchange  # wuh
                    if min == 60:
                        min = 0
                        if it[i] == 23:
                            it[i] = 0
                            leapyear = leap_year(iy[i])
                            if id[i] - leapyear == 365:
                                id[i] = 0
                                iy[i] += 1
                        else:
                            it[i] += 1

                    met_temp[j, 0] = iy[i]  #timestep_now[0]  # iy
                    met_temp[j, 1] = id[i]  #timestep_now[1]  # id
                    met_temp[j, 2] = it[i]  #timestep_now[2]  # it
                    met_temp[j, 3] = min  #timestep_now[3] + min  # imin

                met_new[i * howmanyfives:i * howmanyfives + howmanyfives, :] = met_temp

        return met_new

    def from5minto1hour_v1(self, results):

        # Time columns
        id = results[:, 1]
        it = results[:, 2]
        # imin = results[:, 3]

        suews_1h = np.zeros((results.shape[0], results.shape[1]))
        index = 0
        sumcol = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 50, 51, 52, 53, 62, 63, 64]
        lastcol = [33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 47, 48, 49, 60, 62]
        for i in range(int(id.min()), int(id.max())):
            suews_1h[index, :] = np.mean(results[index * 12: index * 12 + 12, :], axis=0)

            for j in range(0, sumcol.__len__()):
                suews_1h[index, sumcol[j]] = np.sum(results[index * 12: index * 12 + 12, sumcol[j]], axis=0)

            for j in range(0, lastcol.__len__()):
                suews_1h[index, lastcol[j]] = results[index * 12 + 12, lastcol[j]]

            index += 1

        return suews_1h

    def translatemetdata(self, old, ver, inputdata, outputdata, delim):

        # old = old SOLWEIG input file where time is moved (i.e. 12 is between 12 and 13)
        # ver = what version that should be saved (e.g. 2015)
        # inputdata = input text file
        # output data = output text file
        # delim = delimiter of output file

        if old == 1:
            # inputdata = 'M:/SOLWEIG/Inputdata/metdata/gbg_060726.txt'
            met_old = np.loadtxt(inputdata, skiprows=1)
            doy_exist = 0
            dectime_exist = 0
            yyyy = met_old[:, 0]
            mm = met_old[:, 1]
            dd = met_old[:, 2]
            hh = met_old[:, 3]
            Ta = met_old[:, 4].copy()
            RH = met_old[:, 5].copy()
            G = met_old[:, 6].copy()
            D = met_old[:, 7].copy()
            I = met_old[:, 8].copy()
            Ws = met_old[:, 9].copy()

            # Moving one hour
            Ta[1:np.size(Ta)] = Ta[0:np.size(Ta) - 1]
            Ta[0] = Ta[1]
            RH[1:np.size(RH)] = RH[0:np.size(RH) - 1]
            RH[0] = RH[1]
            G[1:np.size(G)] = G[0:np.size(G) - 1]
            G[0] = G[1]
            D[1:np.size(D)] = D[0:np.size(D) - 1]
            D[0] = D[1]
            I[1:np.size(I)] = I[0:np.size(I) - 1]
            I[0] = I[1]
            Ws[1:np.size(Ws)] = Ws[0:np.size(Ws) - 1]
            Ws[0] = Ws[1]
            minute = Ta * 0

        else:
            met_old = np.loadtxt(inputdata, skiprows=1)

            user_input = int(input('Put in manually or translate from v2014 (1 or 0)?: '))
            yyyy_exist = int(input('yyyy exist (1 or 0)?: '))
            if yyyy_exist == 1:
                yyyy_col = int(input('column for yyyy: ')) - 1
            else:
                yy = int(input('Specify year (yyyy): '))

            if user_input == 1:
                doy_exist = int(input('doy exist (1 or 0)?: '))
                if doy_exist == 1:
                    doy_col = int(input('column for doy: ')) - 1
                else:
                    month_col = int(input('column for month: ')) - 1
                    day_col = int(input('column for day of month: ')) - 1

                hh_col = int(input('column for hour: ')) - 1
                dectime_exist = int(input('dectime exist (1 or 0)?: '))
                if dectime_exist == 1:
                    dectime_col = int(input('column for dectime: ')) - 1
                    if ver == 2015:
                        dechour = (met_old[:, dectime_col] - np.floor(met_old[:, dectime_col])) * 24
                        minute = np.round((dechour - np.floor(dechour)) * 60)
                        minute[(minute == 60)] = 0
                else:
                    min_col = int(input('column for min: ')) - 1

                only_mandatory = int(input('Only put in mandatory [Ta,RH,Kdn,pres,Ws] (1 or 0)?: '))

                if only_mandatory == 1:
                    wind_col = int(input('column for Ws: ')) - 1
                    RH_col = int(input('column for RH: ')) - 1
                    Ta_col = int(input('column for Ta: ')) - 1
                    press_exist = int(input('Pressure exist (1 or 0)?: '))
                    if press_exist == 1:
                        press_col = int(input('column for Pressure (kPa): ')) - 1
                    else:
                        press_av = 101.3
                        print 'Pressure set to 101.3 kPa'

                    grad_col = int(input('column for Kdn: ')) - 1
                else:
                    Qstar_col = int(input('column for Q*: ')) - 1
                    Qh_col = int(input('column for Qh: ')) - 1
                    Qe_col = int(input('column for Qe: ')) - 1
                    Qs_col = int(input('column for Qs: ')) - 1
                    Qf_col = int(input('column for Qf: ')) - 1
                    wind_col = int(input('column for Ws: ')) - 1
                    RH_col = int(input('column for RH: ')) - 1
                    Ta_col = int(input('column for Ta: ')) - 1
                    press_exist = int(input('Pressure exist (1 or 0)?: '))
                    if press_exist == 1:
                        press_col = int(input('column for Pressure (kPa): ')) - 1
                    else:
                        press_av = 101.3
                        print 'Pressure set to 101.3 kPa'
                    rain_col = int(input('column for rain: ')) - 1
                    grad_col = int(input('column for Kdn: ')) - 1
                    snow_col = int(input('column for snow: ')) - 1
                    ldown_col = int(input('column for ldown: ')) - 1
                    fcld_col = int(input('column for fcld: ')) - 1
                    wuh_col = int(input('column for wuh: ')) - 1
                    xsmd_col = int(input('column for xsmd: ')) - 1
                    lai_col = int(input('column for lai: ')) - 1
                    drad_col = int(input('column for kdiff: ')) - 1
                    irad_col = int(input('column for kdir: ')) - 1
                    wdir_col = int(input('column for Wdir: ')) - 1

            else:
                doy_col = 1 - 1
                hh_col = 2 - 1
                dectime_col = 3 - 1
                if ver == 2015:
                    dechour = (met_old[:, dectime_col] - np.floor(met_old[:, dectime_col])) * 24
                    minute = np.round((dechour - np.floor(dechour)) * 60)
                    minute[(minute == 60)] = 0

                Qstar_col = 4 - 1
                Qh_col = 5 - 1
                Qe_col = 6 - 1
                Qs_col = 7 - 1
                Qf_col = 8 - 1
                wind_col = 9 - 1
                RH_col = 10 - 1
                Ta_col = 11 - 1
                press_col = 12 - 1
                rain_col = 13 - 1
                grad_col = 14 - 1
                snow_col = 15 - 1
                ldown_col = 16 - 1
                fcld_col = 17 - 1
                wuh_col = 18 - 1
                xsmd_col = 19 - 1
                lai_col = 20 - 1
                drad_col = 20 - 1
                irad_col = 20 - 1
                wdir_col = 20 - 1

            if doy_exist == 1:
                doy = met_old[:, doy_col]
            else:
                mm = met_old[:, month_col]
                dd = met_old[:, day_col]

            if dectime_exist == 0:
                minute = met_old[:, min_col]

            if yyyy_exist == 1:
                yyyy = met_old[:, yyyy_col]
            else:
                yyyy = met_old[:, 1] * 0 + yy

            hh = met_old[:, hh_col]

            if only_mandatory == 0:
                Qstar = met_old[:, Qstar_col]
                Qh = met_old[:, Qh_col]
                Qe = met_old[:, Qe_col]
                Qs = met_old[:, Qs_col]
                Qf = met_old[:, Qf_col]
                rain = met_old[:, rain_col]
                snow = met_old[:, snow_col]
                ldown = met_old[:, ldown_col]
                fcld = met_old[:, fcld_col]
                wuh = met_old[:, wuh_col]
                xsmd = met_old[:, xsmd_col]
                lai= met_old[:, lai_col]
                D = met_old[:, drad_col]
                I = met_old[:, irad_col]
                wdir = met_old[:, wdir_col]

            Ws = met_old[:, wind_col]
            Ta = met_old[:, Ta_col]
            RH = met_old[:, RH_col]

            if press_exist == 1:
                press = met_old[:, press_col]
            else:
                press = met_old[:, hh_col] * 0 + press_av

            G = met_old[:, grad_col]

        if ver == 2015:
            met_new = np.zeros((met_old.shape[0], 24)) - 999
        else:
            met_new = np.zeros((met_old.shape[0], 23)) - 999

        if doy_exist == 0:
            doy = np.zeros((met_old.shape[0]))

        for i in range(0, yyyy.size):
            # day of year and check for leap year
            yy = int(yyyy[i])

            if (yy % 4) == 0:
                if (yy % 100) == 0:
                    if (yy % 400) == 0:
                        leapyear = 1
                    else:
                        leapyear = 0
                else:
                    leapyear = 1
            else:
                leapyear = 0

            if leapyear == 1:
                dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            else:
                dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

            if doy_exist == 0:
                doy[i] = sum(dayspermonth[0:int(mm[i] - 1)]) + dd[i]  # might now work

        if ver == 2015:  # v2015a
            #inputnew = 'M:/SOLWEIG/Inputdata/metdata/2015a/gbg20051011_2015a.txt'
            #f = open(inputnew, 'r')
            #header = f.readline()
            header = '%iy id  it imin   Q*      QH      QE      Qs      Qf    Wind    RH     Td     press   rain ' \
                     '   Kdn    snow    ldown   fcld    wuh     xsmd    lai_hr  Kdiff   Kdir    Wd'

            met_new[:, 0] = yyyy
            met_new[:, 1] = doy
            met_new[:, 2] = hh
            met_new[:, 3] = minute
            met_new[:, 9] = Ws
            met_new[:, 10] = RH
            met_new[:, 11] = Ta
            met_new[:, 14] = G
            met_new[:, 12] = press

            if only_mandatory == 0:
                met_new[:, 4] = Qstar
                met_new[:, 5] = Qh
                met_new[:, 6] = Qe
                met_new[:, 7] = Qs
                met_new[:, 8] = Qf
                met_new[:, 13] = rain
                met_new[:, 15] = snow
                met_new[:, 16] = ldown
                met_new[:, 17] = fcld
                met_new[:, 18] = wuh
                met_new[:, 19] = xsmd
                met_new[:, 20] = lai
                met_new[:, 21] = D
                met_new[:, 22] = I
                met_new[:, 22] = wdir

            # #Save as text files
            numformat = '%3d %2d %3d %2d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f'
            np.savetxt(outputdata, met_new, fmt=numformat, header=header, comments='')

        else:  # v2014a
            inputnew = 'M:/SOLWEIG/Inputdata/metdata/2014a/gbg20051011_2014a_space.txt'
            f = open(inputnew, 'r')
            header = f.readline()

            if dectime_exist == 1:
                dectime = met_old[:, dectime_col]
            else:
                dectime = doy + hh / 24 + minute / (60 * 24)

            met_new[:, 0] = doy
            met_new[:, 1] = hh
            met_new[:, 2] = dectime
            met_new[:, 8] = press
            met_new[:, 9] = RH
            met_new[:, 10] = Ta
            met_new[:, 13] = G
            met_new[:, 20] = D
            met_new[:, 21] = I
            met_new[:, 8] = Ws

            # %Save as text files
            numformat = '%3d %3d %6.5f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f'
            np.savetxt(outputdata, met_new, fmt=numformat, delimiter=delim, header=header, comments='')

