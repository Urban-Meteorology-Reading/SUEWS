__author__ = 'Fredrik Lindberg'

# This class will be used to plot output result from Suews

import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as dt

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


class SuewsPlotting:
    def __init__(self):
        pass

    def make_dectime(self, dataout):
        datenum_yy = np.zeros(dataout.shape[0])
        for i in range(0, dataout.shape[0]): # making date number
            datenum_yy[i] = dt.date2num(dt.datetime.datetime(int(dataout[i, 0]), 1, 1))

        dectime = datenum_yy + dataout[:, 4]

        return dectime

    def plot1hour(self, dataout, datain, dectime):

        plt.figure(figsize=(15, 7), facecolor='white')

        ax1 = plt.subplot(3, 1, 1)
        ax1.plot(dectime, dataout[:, 5], label='$K_{down}$')
        ax1.plot(dectime, dataout[:, 6], label='$K_{up}$')
        ax1.plot(dectime, dataout[:, 7], label='$L_{down}$')
        ax1.plot(dectime, dataout[:, 8], label='$L_{up}$')
        ax1.plot(dectime, dataout[:, 10], label='$Q*$')
        ax1.set_ylim([-100, 1000])
        ax1.set_ylabel(r'$W/m^2$')
        plt.legend()

        ax2 = plt.subplot(3, 1, 2, sharex=ax1)
        ax2.plot(dectime, dataout[:, 13], label='$Q_S$')
        ax2.set_ylabel(r'$W/m^2$')
        ax2.plot(dectime, dataout[:, 14], label='$Q_F$')
        ax2.plot(dectime, dataout[:, 15], label='$Q_H$')
        ax2.plot(dectime, dataout[:, 16], label='$Q_E$')
        ax2.set_ylim([-100, 400])
        plt.legend()

        ax3 = plt.subplot(3, 1, 3, sharex=ax1)
        ax4 = ax3.twinx()
        ax3.plot(dectime, datain[:, 11], 'g-', label='$T_a$')
        ax4.plot(dectime, datain[:, 10], 'b-', label='$RH$')
        ax3.set_xlabel('decimal time')
        ax3.set_ylabel('$T_a (degC)$', color='g')
        ax4.set_ylabel('RH (%)', color='b')
        ax3.set_xlim([min(dectime), max(dectime)])

        plt.show()

    def plotmonthlystatistics(self, dataout, datacol, dectime):

        for i in range(0, datacol.__len__()):
            print i


