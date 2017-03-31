import sys
import os
import math
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar


N = 384000.0
sigma = 1.0

if __name__ == '__main__':

    print('Just GO')




    # probability of noise hits  
    x = np.arange(9500,10000)/10000.
    x = np.arange(500,1000)/1000.

    # number of noise hits for N channels
    n_hits = N*(1.0-x)

    # calculate x of normal distributin for a given probability
    ppf = st.norm.ppf(x)
    

    # CDF of normal distribution
    cdf = st.norm.cdf(x)

    print(x)
    print(ppf)

    fig1, ax1 = plt.subplots(2,1)
    ax1[0].plot(x, ppf)
    ax1[0].set_ylabel('ppf')
    ax1[0].set_xlabel('x')
    ax1[1].plot(x, cdf)
    ax1[1].set_ylabel('cdf')
    ax1[1].set_xlabel('x')

    fig2, ax2 = plt.subplots(3,1)
    fig2.patch.set_facecolor('white')
    ax2[0].set_axis_bgcolor('white')
    ax2[0].plot(ppf/(sigma),x )
    ax2[0].set_xlabel(r'$\sigma$ cut')
    ax2[0].set_ylabel('ppf probability')

    ax2[1].set_axis_bgcolor('white')
    ax2[1].plot(ppf/(sigma),1.0-x )
    ax2[1].set_xlabel(r'$\sigma$ cut')
    ax2[1].set_ylabel('low tail probability')

    ax2[2].set_axis_bgcolor('white')
    ax2[2].semilogy(ppf/(sigma),n_hits )
    ax2[2].set_xlabel(r'$\sigma$ cut')
    ax2[2].set_ylabel('Number of noise hits (384,000 channels)')

    #fig3,ax3 = plt.subplots()
    #ax3.plot(np.arange(-100,100)/10.0/sigma, st.norm.pdf(np.arange(-100,100)/10.0))


    plt.show()
