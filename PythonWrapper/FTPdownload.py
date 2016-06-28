# Funtion:
# Download WATCH/WFDEI data for specified date range and variable.
# Author:
# Lingbo Xue, L.Xue@student.reading.ac.uk
# Ting Sun, ting.sun@reading.ac.uk
# History:
# LX, 13 Jun 2016: initial version
# TS, 14 Jun 2016 : path parsing improved.
# TS, 15 Jun 2016: date/time bugs fixed and logic improved.

from ftplib import FTP
import os
import sys
import urllib
import numpy as np
from datetime import date,datetime
from dateutil.relativedelta import relativedelta


def single_file_download(path, key_word):
    path = os.path.join(path, key_word)  # safely join pathnames
    path = os.path.expanduser(path)  # expand Unix
    if os.path.lexists(path):
        print("File exist! Try again...")
    else:
        f = open(path, 'wb')
        print "%s downloading..." % key_word
        ftp.retrbinary('RETR %s' % key_word, f.write)
        print "%s download succeed!" % key_word


def time_period_files_download(key_word, path):
    ftp.cwd(key_word)
    list = ftp.nlst()

    path = os.path.expanduser(path)  # expand Unix
    os.chdir(path)
    path = os.path.join(path, key_word)  # safely join pathnames
    if not os.path.lexists(path):
        os.mkdir(key_word)

    time1, time2 = 190101, 200112
    time_range = time_period_test(time1, time2)

    if key_word == "Rainf_WFD" or key_word == "Rainf_daily_WFD" or key_word == "Snowf_WFD" or key_word == "Snowf_daily_WFD":
        name_base = key_word + '_CRU'
    else:
        name_base = key_word

    # bug here: such range cannot go over a year: 190113 will be generated.
    for i in time_range:
        file_name = name_base + '_' + str(i) + '.nc'
        # path_temp = path + file_name
        path_temp = os.path.join(path, file_name)
        if os.path.lexists(path_temp):
            print "%s exists!" % file_name
        else:
            f = open(path_temp, 'wb')
            print "%s downloading..." % file_name
            ftp.retrbinary('RETR %s' % file_name, f.write)
            print "%s download succeed!" % file_name


def time_period_test(time1, time2):
    # convert to stirng then parse to date
    time1 = datetime.strptime(str(time1), "%Y%m")
    time2 = datetime.strptime(str(time2), "%Y%m")
    while True:
        print("The time period is from %s to %s" %
              (time1.strftime("%Y%m"), time2.strftime("%Y%m")))
        # get valid date input:
        while True:
            try:
                start_time = input("Please input the start time(YYYYMM): ")
                end_time = input(
                    "Please input the end time(included)(YYYYMM): ")
                # convert to stirng then parse to date
                start_time = datetime.strptime(str(start_time), "%Y%m")
                end_time = datetime.strptime(str(end_time), "%Y%m")
            except:
                print("Sorry, please input a valid time(YYYYMM).")
                # better try again... Return to the start of the loop
                continue

            # check range availability
            if not(time1 <= start_time <= end_time <= time2):
                print("Sorry, please input a valid time range between %s and %s" % (
                    time1.strftime("%Y%m"), time2.strftime("%Y%m")))
                continue
            else:
                # valid range received
                # we're ready to exit the loop.
                break

        # add valid months to a list
        print("the time range between %s and %s" %
              (start_time.strftime("%Y%m"), end_time.strftime("%Y%m")))
        range_time = [start_time.strftime("%Y%m")]
        while start_time < end_time:
            start_time = start_time + relativedelta(months=+1)
            range_time.append(start_time.strftime("%Y%m"))
        return range_time

# login
ftp = FTP("ftp.iiasa.ac.at")
ftp.login("rfdata", "forceDATA")
print ftp.getwelcome()

# download
while True:
    ftp.cwd("/WATCH_Forcing_Data")
    print ftp.dir()

    while True:
        try:
            key_word = raw_input(
                "Please choose the data you want to download from the list above: ")
            ftp.nlst(key_word)
        except:
            print "No such file or directory. Try again..."
        else:
            break

    while True:
        path = raw_input("Please input the download path: ")
        path = os.path.realpath(os.path.expanduser(path))
        print path
        if os.path.lexists(path):
            break
        else:
            print "No such directory. Try again..."
        # /Users/xuelingbo/Desktop/test/
        print key_word

    if '.' in key_word:
        single_file_download(path, key_word)
    elif key_word == "Rainf_WFD" or key_word == "Snowf_WFD":
        ftp.cwd(key_word)
        list = ftp.nlst()
        if input("Download %s [1]or files from 190101 to 200112[2]? 1/2" % list[0]) == 1:
            os.chdir(path)
            path = path + key_word + '/'
            if not os.path.lexists(path):
                os.mkdir(key_word)
            single_file_download(path, list[0])
        else:
            ftp.cwd("/WATCH_Forcing_Data")
            time_period_files_download(key_word, path)
    elif key_word == "WFDEI":
        ftp.cwd(key_word)
        print ftp.dir()

        os.chdir(path)
        # path += key_word + '/'
        path = os.path.join(path, key_word)  # safely join pathnames
        path = os.path.realpath(os.path.expanduser(path))
        if not os.path.lexists(path):
            os.mkdir(key_word)

        while True:
            try:
                key_word_2 = raw_input(
                    "Please choose the data you want to download from the list above: ")
                ftp.nlst(key_word_2)
            except:
                print "No such file or directory. Try again..."
            else:
                break

        if '.' in key_word_2:
            single_file_download(path, key_word_2)
        else:
            ftp.cwd(key_word_2)

            os.chdir(path)
            # path += key_word_2 + '/'
            path = os.path.join(path, key_word_2)  # safely join pathnames
            path = os.path.realpath(os.path.expanduser(path))

            print path
            print os.path.lexists(path)
            if not os.path.lexists(path):
                os.mkdir(key_word_2)
            # os.chdir(path)

            if key_word_2[-4:] == "GPCC":
                time1, time2 = 197901, 201312
            else:
                time1, time2 = 197901, 201412
            time_range = time_period_test(time1, time2)

            name_base = key_word_2

            for i in time_range:
                file_name = name_base + '_' + str(i) + '.nc.gz'
                # path_temp = path + file_name
                # safely join pathnames
                path_temp = os.path.join(path, file_name)

                if os.path.lexists(path_temp):
                    print "%s exists!" % file_name
                else:
                    f = open(path_temp, 'wb')
                    print "%s downloading..." % file_name
                    ftp.retrbinary('RETR %s' % file_name, f.write)
                    print "%s download succeed!" % file_name
    else:
        time_period_files_download(key_word, path)

    t = raw_input('Do you want to quit? Y/N')
    if t == 'Y' or t == 'y':
        ftp.quit()
        break
