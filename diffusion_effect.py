
import sys
import os
import math
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import diffusion



def get_args():
    parser = argparse.ArgumentParser('Signal analysis')
    parser.add_argument('file', nargs=1, help='Data file.')
    args = parser.parse_args()
    print(args)
    return args


def gaus(x,mu,sig):
    return 1.0/np.sqrt(2*np.pi*np.power(sig,2))*1.0/np.exp(np.power(x-mu,2)/(2*np.power(sig,2)))


if __name__ == '__main__':
    print('Just GO')

    args = get_args()


    
    Ts = 0.1e-6 # sampling interval in seconds
    Fs = 1.0/Ts # sampling rate in Hz
    
    # find histogram of time signal    
    a = np.loadtxt(args.file[0])
    y = a[:,1]
    t = a[:,0]*1e-6 # turn into seconds


    #y = y[n/2:]
    #t = t[len(t)/2:]

    n = len(y) # length of signal
    k = np.arange(n)
    T = n/Fs # length in time of signal
    frq = k/T # two side frequency range
    frq = frq[range(n/2)] # one side frequency range

    Y = np.fft.fft(y) # /n FFT calculation
    Yp = Y[range(1,(n-1)/2+1)] # one sided
    Ym = Y[range((n+1)/2,n)] # one sided

    

    


    for d in [0.001,100.,200.,300.,360.]:

        saveFilename = os.path.splitext(os.path.basename(args.file[0]))[0] + '_d{0:.3f}cm'.format(d)

        td = diffusion.diffusion_time(d)

        #make diffusion gaussian
        t_g = np.arange(0,100)*Ts - 100*Ts/2
        g = gaus(t_g,0,td)
        g = g/np.sum(g)

        # plit gaussian
        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)
        ax0.plot(t_g*1e6,g) # time signal
        ax0.set_title('Diffusion time')
        ax0.set_xlabel(r'Time ($\mu$s)')
        ax0.set_ylabel(r'Arbitrary units')
        ax0.text(0.01, 0.9, r'drift distance {0:f}cm'.format(d), transform=ax0.transAxes, color='green', fontsize=10)
        fig0.savefig(saveFilename + '_diffusion.png',bbox_inches='tight')

        # convolve
        yd = np.convolve(y,g)
        nd = len(yd)

        # convolved time axis is longer
        #t_ext = np.concatenate(t,
        t_ext= np.arange(nd-n)*Ts + Ts + t[-1]
        print('t ' + str(t))
        print('t_ext ' + str(t_ext))
        print(t_ext)
        print(t.shape)
        print(t_ext.shape)
        print(str(n))
        print(str(nd))
        print(str(nd-n))
        print(yd)
        print(td)
        print(str(t[len(t)/2]))
        t_ext2 = np.concatenate((t,t_ext))
        
        
        # calculate FFT
        Yd = np.fft.fft(yd) # /n FFT calculation
        Ydp = Yd[range(nd/2)] # one sided
        Ydm = Yd[range(nd/2,nd)] # one sided        
        T = nd/Fs # length in time of signal
        frq_d = np.arange(nd)/T # two side frequency range
        frq_d = frq_d[range(nd/2)] # one side frequency range

        # plot time signal and frequency spectrum of input signal
        fig1, ax1 = plt.subplots(2,1)
        ax1[0].plot(t,y) # time signal
        ax1[0].set_label('Raw input signal')
        ax1[0].set_xlabel(r'Time ($\mu$s)')
        ax1[0].set_ylabel(r'Current response ($\mu$A)')
        ax1[1].semilogx(frq,np.abs(Yp)) #color='red') # amplitude freq.
        ax1[1].set_xlabel('Frequency spectrum [Hz]')
        ax1[1].set_ylabel('Amplitude')
        ax1[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax1[0].transAxes, color='green', fontsize=10)
        #ax1[0].text(0.01, 0.8, r'd={0:f}cm'.format(d), transform=ax1[0].transAxes, color='green', fontsize=10)
        fig1.savefig(saveFilename + '_input_signal.png',bbox_inches='tight')

        # plot time signal and frequency spectrum of diffused input signal

        fig2, ax2 = plt.subplots(2,1)
        ax2[0].plot(t_ext2, yd,color='r') # time signal
        ax2[0].set_label('Diffused input signal')
        ax2[0].set_xlabel(r'Time ($\mu$s)')
        ax2[0].set_ylabel(r'Current response ($\mu$A)')
        ax2[1].semilogx(frq_d,np.abs(Ydp),color='r') #color='red') # amplitude freq.
        ax2[1].set_xlabel('Frequency spectrum [Hz]')
        ax2[1].set_ylabel('Amplitude')
        ax2[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax2[0].transAxes, color='green', fontsize=10)
        ax2[0].text(0.01, 0.8, r'd={0:f}cm'.format(d), transform=ax2[0].transAxes, color='green', fontsize=10)
        fig2.savefig(saveFilename + '_diffused_input_signal.png',bbox_inches='tight')


        #plt.show()
