
import sys
import os
import math
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar



def get_args():
    parser = argparse.ArgumentParser('Signal analysis')
    parser.add_argument('--tau','-t', type=float, required=True, help='Filter time constant in micro seconds.')
    parser.add_argument('--freq', type=float, default=1.0, help='ADC frequency in Msps.')
    parser.add_argument('file', nargs=1, help='Data file.')
    parser.add_argument('--batch', action='store_true')
    args = parser.parse_args()
    print(args)
    return args
    


if __name__ == '__main__':
    print('Just GO')

    args = get_args()
    
    tau = args.tau*1e-6 

    saveFilename = os.path.splitext(os.path.basename(args.file[0]))[0] + '_{0:.3f}usec_{1:.3f}Msps'.format(args.tau,args.freq)

    Ts = 0.1e-6 # sampling interval in seconds
    Fs = 1.0/Ts # sampling rate in Hz
    
    # find histogram of time signal    
    a = np.loadtxt(args.file[0])
    y = a[:,1]
    t = a[:,0]*1e-6 # turn into seconds


    #y = y[len(y)/2:]
    #t = t[len(t)/2:]


    n = len(y) # length of signal
    k = np.arange(n)
    T = n/Fs # length in time of signal
    frq = k/T # two side frequency range
    frq = frq[range(n/2)] # one side frequency range

    # calculate FFT
    Y = np.fft.fft(y) # /n FFT calculation
    Yp = Y[range(1,(n-1)/2+1)] # one sided
    Ym = Y[range((n+1)/2,n)] # one sided

    # plot time signal and frequency spectrum of input signal
    fig1, ax1 = plt.subplots(2,1)
    ax1[0].plot(t,y) # time signal
    ax1[0].set_xlabel('Time (microsec)')
    ax1[0].set_ylabel('Current response (microA)')
    ax1[1].plot(frq,np.abs(Yp),'r') #color='red') # amplitude freq.
    ax1[1].plot(frq,np.flipud(np.abs(Ym)),'c') #color='red') # amplitude freq.
    ax1[1].set_xlabel('Frequency spectrum [Hz]')
    ax1[1].set_ylabel('Amplitude')
    ax1[0].set_title('Input signal')
    ax1[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax1[0].transAxes, color='green', fontsize=10)
    fig1.savefig(saveFilename + '_input_signal.png',bbox_inches='tight')
    

    # 1st order low-pass filter (RC)
    #tau = 5.0e-6 #time constant
    fc = 1/(2*math.pi*tau) # cut-off frequency
    t_ir = np.arange(0,100)*Ts # time axis    
    y_ir = np.exp(-t_ir/tau) # impulse response
    #y_ir = np.exp(-t_ir/tau) # impulse response
    n_ir = len(y_ir) # number of samples
    Ts_ir = tau # sampling time
    Fs_ir = 1/Ts_ir # sampling rate
    T_ir = n_ir/Fs_ir # total length
    frq_ir = np.arange(n_ir)/T_ir # frequency range    
    K_ir = 2*math.pi*frq_ir*tau    # helper constant
    gc_ir = 1/np.sqrt(1 + np.square(K_ir))     # direct transfer function for filter
    # use same freq range as signal fft to apply
    K = 2*math.pi*frq*tau    # helper constant
    gc = 1/np.sqrt(1 + np.square(K))     # direct transfer function for filter

    # plot impulse response and transfer function
    fig2, ax2 = plt.subplots(2,1)
    ax2[0].plot(t_ir,y_ir) 
    ax2[0].set_xlabel('Time (microsec)')
    ax2[0].set_ylabel('Arb. Units.')
    ax2[0].set_title('Filter Impulse Response')
    
    ax2[1].plot(frq_ir,gc_ir) 
    #ax2[1].semilogx(frq,gc) # freq
    #ax2[1].semilogx(frq,np.abs(Y_ir)) # freq
    ax2[1].plot([fc,fc],[0,1], color='k', linestyle='-', linewidth=2) 
    ax2[1].plot([fc/10,fc*10],[1/math.sqrt(2),1/math.sqrt(2)], color='k', linestyle='-', linewidth=2) 
    ax2[1].set_xlabel('Frequency spectrum [Hz]')
    ax2[1].set_ylabel('Amplitude')
    ax2[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax2[0].transAxes, color='green', fontsize=10)
    ax2[0].text(0.01, 0.8, r'$\tau$={0:.2f}$\mu$s'.format(tau*1e6), transform=ax2[0].transAxes, color='green', fontsize=10)



    fig2.savefig(saveFilename + '_filter.png',bbox_inches='tight')



    
    # convolve impulse respone with input signal
    y_c = np.convolve(y_ir,y)
    y_c = y_c*np.sum(y)/np.sum(y_c) # normalize to input signal
    Y_C = np.fft.fft(y_c) # do FFT to get frequency spectrum
    Yp_C = Y_C[range(len(Y_C)/2)] # one sided range
    c = (np.arange(len(y_c)-n)*Ts + Ts + t[-1]) # helper constant

    # create time axis of convoluted signal (longer by length of convolution)
    t_y_c = np.concatenate((t,c))
    T_c  = len(t_y_c)*Ts # length of convoluted signal
    frq_c = np.arange(len(t_y_c))/T_c # frequency bins

    # apply transfer function to frequency domain signal representation
    # not sure how to handle neg frequencies here
    Yp_f = Yp*gc 
    Ym_f = np.flipud(np.flipud(Ym)*gc) 
    Y_f = np.concatenate((Yp_f,Ym_f))
    Y_f = np.insert(Y_f,0,Y[0])


    # transform back the signal to time domain
    iY = np.fft.ifft(Y_f)


    # plot comparison of input and filtered by applying transfer function directly to freq. spectrum of signal
    fig4, ax4 = plt.subplots(2,1)
    ax4[0].plot(t,y, color='black') 
    ax4[0].plot(t,iY, color='red') 
    ax4[0].set_xlabel('Time')
    ax4[0].set_ylabel('Arb. Units.')
    ax4[0].set_title('Transfer function filtered response')

    ax4[1].plot(frq,np.abs(Yp), color='black') 
    iYp = iY[range(1,(n-1)/2+1)] # one sided
    ax4[1].plot(frq,np.abs(iYp), color='red') #color='red') # amplitude freq.
    ax4[1].set_title('Transfer function filtered response')
    ax4[1].set_xlabel('Frequency [Hz]')
    ax4[1].set_ylabel('Arb. Units.')

    ax4[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax4[0].transAxes, color='green', fontsize=10)
    ax4[0].text(0.01, 0.8, r'$\tau$={0:.2f}$\mu$s'.format(tau*1e6), transform=ax4[0].transAxes, color='green', fontsize=10)


    fig4.savefig(saveFilename + '_tffiltered_signal.png',bbox_inches='tight')


    # plot comparison of input and filtered using convolution
    fig3, ax3 = plt.subplots(2,1)
    ax3[0].plot(t,y, color='black') 
    ax3[0].plot(t_y_c,y_c, color='red') 
    ax3[0].set_xlabel('Time')
    ax3[0].set_ylabel('Arb. Units.')
    ax3[0].set_title('Filtered response')

    ax3[1].plot(frq,np.abs(Yp), color='black') 
    ax3[1].plot(frq_c[range(len(t_y_c)/2)],np.abs(Yp_C), color='red') 
    ax3[1].set_xlabel('Frequency [Hz]')
    ax3[1].set_ylabel('Arb. Units.')

    ax3[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax3[0].transAxes, color='green', fontsize=10)
    ax3[0].text(0.01, 0.8, r'$\tau$={0:.2f}$\mu$s'.format(tau*1e6), transform=ax3[0].transAxes, color='green', fontsize=10)
                           
    fig3.savefig(saveFilename + '_filtered_signal.png',bbox_inches='tight')



    # sample filtered signal
    print('Sampling f ' + str(args.freq) + ' MHz')
    print('T_c ' + str(T_c) + ' sec total length of conv signal')
    print('n_c ' + str(len(t_y_c)) )
    print('Fs ' + str(Fs) )
    Fs_frac = args.freq/(Fs*1e-6)
    print('Fs_frac ' + str(Fs_frac) )
    adc_samples = np.arange(len(y_c))
    print('adc_samples ' + str(adc_samples))
    adc_samples = adc_samples[ adc_samples*Fs_frac%1==0 ]
    print('adc_samples ' + str(adc_samples))
    y_c_samples = y_c[ adc_samples ]
    t_y_c_samples = t_y_c[ adc_samples ]
    print('y_c ' + str(y_c))
    print('y_c_samples ' + str(y_c_samples))

    # do a FFT on sampled signal
    Y_c_samples = np.fft.fft(y_c_samples)
    Yp_c_samples = Y_c_samples[range(len(Y_c_samples)/2)]
    Ym_c_samples = Y_c_samples[range(len(Y_c_samples)/2,len(Y_c_samples))]
    Fs_samples = Fs*Fs_frac
    T_c_samples = len(y_c_samples)/Fs_samples # length in time of signal
    frq_samples = np.arange(len(y_c_samples))/T_c_samples # freq axis

    # deconvolute with filter transfer function
    K_samples = 2*math.pi*frq_samples*tau     # helper constant, use same freq range as signal fft to apply
    gc_samples = 1/np.sqrt(1 + np.square(K_samples))     # direct transfer function for filter
    Yp_c_samples_deconv = Yp_c_samples/gc_samples[range(len(y_c_samples)/2)] # deconvolve
    Ym_c_samples_deconv = np.flipud(np.flipud(Ym_c_samples)/gc_samples[range(len(y_c_samples)/2)])
    Y_c_samples_deconv = np.concatenate((Yp_c_samples_deconv,Ym_c_samples_deconv))
    




    # inverse FFT to time domain
    
    y_c_samples_deconv = np.fft.ifft(Y_c_samples_deconv)
    y_c_samples_deconv = y_c_samples_deconv*np.sum(y_c_samples)/np.sum(y_c_samples_deconv) # normalize



    print(' len(Y_c_samples) ' + str( len(Y_c_samples)))
    print(' Fs ' + str( Fs ))
    print(' Fs_frac ' + str( Fs_frac) )
    print(' Fs_samples ' + str( Fs_samples) )
    print(' frq ' + str( frq_samples) )
    print(' frq_c ' + str( frq_samples) )
    print(' frq_samples ' + str( frq_samples) )
    print(' y_c_samples.shape ' + str( y_c_samples.shape))
    print(' Y_c_samples.shape ' + str( Y_c_samples.shape))
    print(' Yp_c_samples.shape ' + str( Yp_c_samples.shape) )
    print(' Ym_c_samples.shape ' + str( Ym_c_samples.shape) )
    print(' Yp_c_samples_deconv.shape ' + str( Yp_c_samples_deconv.shape) )
    print(' Ym_c_samples_deconv.shape ' + str( Ym_c_samples_deconv.shape) )
    print(' Y_c_samples_deconv.shape ' + str( Y_c_samples_deconv.shape) )
    print(' y_c_samples_deconv.shape ' + str( y_c_samples_deconv.shape) )


    
    # plot comparison of input and filtered using convolution
    fig5, ax5 = plt.subplots(2,1)
    ax5[0].plot(t,y, color='black') 
    ax5[0].plot(t_y_c,y_c, color='red') 
    ax5[0].plot(t_y_c_samples,y_c_samples, color='blue', linestyle='-', marker='o') 
    ax5[0].plot(t_y_c_samples,np.abs(y_c_samples_deconv), color='green')
    ax5[0].set_xlabel('Time')
    ax5[0].set_ylabel('Arb. Units.')
    ax5[0].set_title('ADC Response')

    #ax5[1].plot(range(len(Y_c_samples)),np.abs(Y_c_samples), color='blue') 
    ax5[1].plot(frq_samples[range(len(y_c_samples)/2)],np.abs(Yp_c_samples), color='blue') 
    ax5[1].plot(frq_samples[range(len(y_c_samples)/2)],np.abs(Yp_c_samples_deconv), color='green') 
    #ax5[1].plot(frq_samples[range(len(y_c_samples)/2)],gc_samples[range(len(y_c_samples)/2)], color='red') 
    #ax5[1].plot(frq_samples[range(len(y_c_samples)/2)],np.abs(Yp_c_samples_deconv), color='black') 
    ax5[1].set_xlabel('Frequency [Hz]')
    ax5[1].set_ylabel('Arb. Units.')

    
    ax5[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax5[0].transAxes, color='green', fontsize=10)
    ax5[0].text(0.01, 0.8, r'$\tau$={0:.2f}$\mu$s'.format(tau*1e6), transform=ax5[0].transAxes, color='green', fontsize=10)
    ax5[0].text(0.01, 0.7, 'ADC {0:.2f}Msps'.format(args.freq), transform=ax5[0].transAxes, color='green', fontsize=10)
                           
    fig5.savefig(saveFilename + '_filtered_signal_sampled.png',bbox_inches='tight')
    

    # plot comparison of input and recovered ADC signal
    fig6 = plt.figure(6)
    ax6 = plt.subplot(111)
    ax6.plot(t,y, color='black') 
    ax6.plot(t_y_c_samples,np.abs(y_c_samples_deconv), color='green', marker='o')
    ax6.set_xlabel('Time')
    ax6.set_ylabel('Current response (microA)')
    ax6.set_title('ADC Deconvoluted Signal')
    ax6.text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax1[0].transAxes, color='green', fontsize=10)
    ax6.text(0.01, 0.8, r'$\tau$={0:.2f}$\mu$s'.format(tau*1e6), transform=ax6.transAxes, color='green', fontsize=10)
    ax6.text(0.01, 0.7, 'ADC {0:.2f}Msps'.format(args.freq), transform=ax6.transAxes, color='green', fontsize=10)
    fig6.savefig(saveFilename + '_adc_deconv_signal.png',bbox_inches='tight')






    if not args.batch:
        plt.show()

    print('Done processing' + args.file[0])

