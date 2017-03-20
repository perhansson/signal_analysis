#!/bin/python

import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar



_E  = 500.0#V/cm
_eL = 2e-2 #eV
_dL = 7.2 #cm^2/s
_v = 0.2*1e6 #cm/s

def diffusion_time(d):
    """Diffusion time for given drift distance for DUNE LArTPC."""
    return np.sqrt(2*d*_dL/(_v*_v*_v))

if __name__ == '__main__':

    print('Just GO')

    print('eL = ' + str(_eL))
    print('dL = ' + str(_dL))

    # drift length
    d = np.arange(0,1000)*0.360 #cm
    #t = d/_v #sec

    # diffusion time
    sigmaL = diffusion_time(d)#np.sqrt(2*d*_dL/(_v*_v*_v))
    sigmaL2 = np.sqrt(2*_eL*d/_E)*1/_v

    # translate to distance
    sigmaD = sigmaL*_v
    sigmaD2 = sigmaL2*_v

    fig1, ax1 = plt.subplots(2,1)

    ax1[0].plot(d,sigmaL*1e6,color='black')
    ax1[0].plot(d,sigmaL2*1e6,color='red')
    ax1[0].set_xlabel('Drift distance (cm)')
    ax1[0].set_ylabel('Diffusion time RMS (us)')

    ax1[0].text(0.01, 0.9, 'arxiv 1508.07059', transform=ax1[0].transAxes, color='green')

    ax1[1].plot(d,sigmaD,color='black')
    ax1[1].plot(d,sigmaD2,color='red')
    ax1[1].set_xlabel('Drift distance (cm)')
    ax1[1].set_ylabel('Spat. res. from diffusion time (cm)')


    plt.show()


            
