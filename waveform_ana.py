
import sys
import os
import math
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import diffusion
import diffusion_effect



def get_args():
    parser = argparse.ArgumentParser('Signal analysis')
    parser.add_argument('--tau','-t', type=float, required=True, help='Filter time constant in micro seconds.')
    parser.add_argument('--freq', type=float, default=1.0, help='ADC frequency in Msps.')
    parser.add_argument('file', nargs=1, help='Data file.')
    parser.add_argument('--batch', action='store_true')
    parser.add_argument('--distance','-d', default=0.0, type=float, help='Drift distance in cm.')
    parser.add_argument('--diffusion', action='store_true', help='Add longitudinal diffusion to waveform.')
    parser.add_argument('--filter', default='rc', help='Filter type')
    parser.add_argument('--SN', default=9.0, type=float, help='Signal to noise ratio')
    args = parser.parse_args()
    print(args)
    return args
    

def gaus(*args):
    return 1.0/np.exp( np.power(args[0],2)/(2*np.power(args[1],2)) )

def rc(*args):
    return np.exp(-args[0]/args[1])

def rc_tf(*args):
    return 1/np.sqrt(1 + np.square(2.0*np.pi*args[0]*args[1]))

def gaus_tf(*args):
    return 1.0/np.exp( np.power(args[0],2) / (2*np.power(args[1],2)) )



if __name__ == '__main__':
    print('Just GO')

    args = get_args()
    
    tau = args.tau*1e-6 

        
    saveFilename = os.path.splitext(os.path.basename(args.file[0]))[0] + '_{0:.3f}usec_{1:.3f}Msps_diff{2:.3f}cm_SN{3:.1f}'.format(args.tau,args.freq, args.distance, args.SN)

    Ts = 0.1e-6 # sampling interval in seconds
    Fs = 1.0/Ts # sampling rate in Hz
    
    # find histogram of time signal    
    a = np.loadtxt(args.file[0])
    y = a[:,1]
    t = a[:,0]*1e-6 # turn into seconds
    n = len(y) # length of signal
    k = np.arange(n)
    T = n/Fs # length in time of signal
    frq = k/T # two side frequency range
    frq = frq[range(n/2)] # one side frequency range

    # calculate FFT of input signal
    Y = np.fft.fft(y) # /n FFT calculation
    if n%2==0:
        max_pos = n/2
        max_neg = n/2
    else:
        max_pos = (n-1)/2
        max_neg = (n+1)/2
    Yp = Y[range(1,max_pos+1)] # one sided pos freq range
    Ym = Y[range(max_neg,n)] # one sided neg freq range

    

    # plot time signal and frequency spectrum of input signal
    fig1, ax1 = plt.subplots(2,1)
    ax1[0].plot(t*1e6,y) # time signal
    ax1[0].set_xlabel('Time (microsec)')
    ax1[0].set_ylabel('Current response (microA)')
    ax1[1].semilogx(frq,np.abs(Yp), color='red') # amplitude freq.
    #ax1[1].plot(frq,np.flipud(np.abs(Ym)), color='pink') # amplitude freq.
    ax1[1].set_xlabel('Frequency spectrum [Hz]')
    ax1[1].set_ylabel('Amplitude')
    ax1[0].set_title('Input signal')
    ax1[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax1[0].transAxes, color='green', fontsize=10)
    fig1.savefig(saveFilename + '_input_signal.png',bbox_inches='tight')


    # add diffusion
    if args.diffusion and args.distance>0.:

        # find longitudinal diffusion time
        time_diffusion = diffusion.diffusion_time(args.distance)
        print('time_diffusion ' + str(time_diffusion))
        
        # make diffusion Gaussian
        #t_d = np.arange(0,100)*Ts - 100*Ts/2
        t_d = np.arange(0,100)*time_diffusion*10 - 100*time_diffusion*10/2
        t_d/=100.0
        y_g = diffusion_effect.gaus(t_d,0,time_diffusion)
        y_g = y_g/np.sum(y_g)

        # plot diffusion gaussian
        fig11, ax11 = plt.subplots(2,1)
        ax11[0].plot(t_d*1e6,y_g) # time signal
        ax11[0].set_title('Diffusion time')
        ax11[0].set_xlabel(r'Time ($\mu$s)')
        ax11[0].set_ylabel(r'Arbitrary units')
        ax11[0].text(0.01, 0.9, r'drift distance {0:f}cm'.format(args.distance), transform=ax11[0].transAxes, color='green', fontsize=10)
        ax11[0].text(0.01, 0.8, r'diff. time $\sigma$={0:.2f}$\mu$s'.format(time_diffusion*1e6), transform=ax11[0].transAxes, color='green', fontsize=10)

        # convolve input signal with diffusion
        y_d = np.convolve(y, y_g)

        # create a new time axis to match the new length
        c = (np.arange(len(y_d)-n)*Ts + Ts + t[-1]) # helper constant
        t_y_d = np.concatenate((t,c))
        T_d  = len(t_y_d)*Ts # length of convoluted signal
        frq_d = np.arange(len(t_y_d))/T_d # frequency bins

        print('y ' + str(y))
        print('sum(y) = ' + str(np.sum(y)))
        print('len(y) ' + str(len(y)))
        print('y_g ' + str(y_g))
        print('len(y_g) ' + str(len(y_g)))
        print('y_d ' + str(y_d))
        print('len(y_d) ' + str(len(y_d)))
        print('t_y_d ' + str(t_y_d))
        print('len(t_y_d) ' + str(len(t_y_d)))

        ax11[1].plot(t_y_d*1e6, y_d) # time signal
        ax11[1].plot(t*1e6, y) # time signal
        ax11[1].set_title('Diffused signal')
        ax11[1].set_xlabel(r'Time ($\mu$s)')
        ax11[1].set_ylabel(r'Arbitrary units')
        ax11[1].text(0.01, 0.9, r'drift distance {0:f}cm'.format(args.distance), transform=ax11[1].transAxes, color='green', fontsize=10)

        fig11.savefig(saveFilename + '_diffusion.png',bbox_inches='tight')

        print('t_d ' + str(t_d))
        print('y_d ' + str(y_d))
        
        # use the new time axis going forward
        t = t_y_d
        frq = frq_d
        y = y_d
        n = len(y_d) # length of signal
        k = np.arange(n)
        T = n/Fs # length in time of signal
        frq = k/T # two side frequency range
        frq = frq[range(n/2)] # one side frequency range

        # calculate FFT of input signal
        Y = np.fft.fft(y) # /n FFT calculation
        
        if n%2==0:
            max_pos = n/2
            max_neg = n/2
        else:
            max_pos = (n-1)/2
            max_neg = (n+1)/2
        Yp = Y[range(1,max_pos+1)] # one sided pos freq range
        Ym = Y[range(max_neg,n)] # one sided neg freq range

        print('n ' + str(n))
        print('len(Y) ' + str(len(Y)))
        print('len(Yp) ' + str(len(Yp)))
        print('len(Ym) ' + str(len(Ym)))
        print('frq ' + str(len(frq)))
        print('frq ' + str(frq))

    

    # plot time signal and frequency spectrum of diffused signal
    fig111, ax111 = plt.subplots(2,1)
    ax111[0].plot(t*1e6,y) # time signal
    ax111[0].set_xlabel('Time (microsec)')
    ax111[0].set_ylabel('Current response (microA)')
    ax111[1].semilogx(frq,np.abs(Yp), color='red') # amplitude freq.
    #ax111[1].semilogy(frq,np.flipud(np.abs(Ym)), color='pink') # amplitude freq.
    ax111[1].set_xlabel('Frequency spectrum [Hz]')
    ax111[1].set_ylabel('Amplitude')
    ax111[0].set_title('Input signal')
    ax111[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax111[0].transAxes, color='green', fontsize=10)
    ax111[0].text(0.01, 0.8, r'drift distance {0:f}cm'.format(args.distance), transform=ax111[0].transAxes, color='green', fontsize=10)
    fig111.savefig(saveFilename + '_input_signal_diffused.png',bbox_inches='tight')

    print('sum(y) = ' + str(np.sum(y)))




    ### Noise

    if 1==1:

        # add white noise: flat in frequency space

        amplitude = np.max(y)
        print('amplitude ' + str(amplitude))
        noise_sigma = amplitude/args.SN
        print('noise_sigma ' + str(noise_sigma))    
        y_noise = np.random.normal(0,noise_sigma,len(y))
        print('y_noise ' + str(y_noise))    

        # plot time signal and frequency spectrum of diffused signal
        fig1111, ax1111 = plt.subplots(2,1)
        ax1111[0].plot(t*1e6,y, label='Signal') 
        ax1111[0].plot(t*1e6,y_noise, label='Noise', color='red')
        ax1111[0].set_xlabel('Time (microsec)')
        ax1111[0].set_ylabel('Current response (microA)')
        #ax1111[1].semilogx(frq,np.abs(Yp), color='red') # amplitude freq.
        #ax111[1].semilogy(frq,np.flipud(np.abs(Ym)), color='pink') # amplitude freq.
        Y_noise = np.fft.fft(y_noise)
        ax1111[1].semilogx(frq,np.abs(Yp),label='Signal') # amplitude freq.
        ax1111[1].semilogx(frq,np.abs(Y_noise)[1:len(y_noise)/2+1], label='Noise') 

        # add noise to signal
        y = y + y_noise
        Y = np.fft.fft(y)
        ax1111[0].plot(t*1e6,y, label='Signal+Noise', color='green')
        ax1111[1].semilogx(frq,np.abs(Y)[1:len(y)/2+1], label='Signal+Noise') 
        legend = ax1111[0].legend(loc='upper left',shadow=True)
        legend = ax1111[1].legend(loc='upper right',shadow=True)
        fig111.savefig(saveFilename + '_input_signal_noise.png',bbox_inches='tight')


    #### Filter
    if args.filter == 'gaus':
        imp = gaus
        tf = gaus_tf
        xmax = 5*tau
        nbins = int(2*xmax/Ts)
        #t_ir = np.arange(0,100)*Ts - 50*Ts # time axis    
        t_ir = np.arange(0,nbins)*Ts - nbins*Ts/2 # time axis    
    elif args.filter == 'rc':
        imp = rc
        tf = rc_tf
        xmax = 5*tau
        nbins = int(xmax/Ts)
        t_ir = np.arange(0,nbins)*Ts # time axis    
        #t_ir = np.arange(0,100)*Ts # time axis    
    else:
        print('invalid filter: ' + args.filter)

    fc = 1/(2*math.pi*tau) # cut-off frequency
    #y_ir = np.exp(-t_ir/tau) # impulse response
    print('Ts ' + str(Ts))
    print('tau ' + str(tau))
    print('t_ir)')
    print(t_ir)
    if args.filter == 'gaus':
        y_ir = imp(t_ir,tau) # impulse response
    elif args.filter == 'rc':
        y_ir = imp(t_ir,tau) # impulse response        
    print('y_ir)')
    print(y_ir)    
    n_ir = len(y_ir) # number of samples
    Ts_ir = Ts # sampling time
    Fs_ir = 1/Ts_ir # sampling rate
    T_ir = n_ir/Fs_ir # total length
    frq_ir = np.arange(n_ir)/T_ir # frequency range    
    print('frq_r')
    print(frq_ir)
    if args.filter == 'gaus':
        gc_ir = tf(frq_ir,1.0/(2*np.pi*tau)) # direct transfer function for filter
        gc = tf(frq,1.0/(2*np.pi*tau))
    elif args.filter == 'rc':
        gc_ir = tf(frq_ir,tau) # direct transfer function for filter        
        gc = tf(frq,tau) # direct transfer function for filter        
    #gc_ir = 1/np.sqrt(1 + np.square(K_ir))     # direct transfer function for filter
    # use same freq range as signal fft to apply
    print('frq')
    print(frq)
    #gc = 1/np.sqrt(1 + np.square(K))     # direct transfer function for filter

    # plot impulse response and transfer function
    fig2, ax2 = plt.subplots(2,1)
    ax2[0].plot(t_ir,y_ir) 
    ax2[0].set_xlabel('Time (microsec)')
    ax2[0].set_ylabel('Arb. Units.')
    ax2[0].set_title('Filter Impulse Response')
    print('gc_ir')
    print(gc_ir)
    
    ax2[1].semilogx(frq_ir,gc_ir) 
    #ax2[1].semilogx(frq,gc) # freq
    #ax2[1].semilogx(frq,np.abs(Y_ir)) # freq
    #ax2[1].plot([fc,fc],[0,1], color='k', linestyle='-', linewidth=2) 
    #ax2[1].plot([fc/10,fc*10],[1/math.sqrt(2),1/math.sqrt(2)], color='k', linestyle='-', linewidth=2) 
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
    if n%2==0:
        Y_f = np.concatenate((Yp_f,Ym_f[1:len(Ym_f)])) # max frequency at the same bin ?
    else:
        Y_f = np.concatenate((Yp_f,Ym_f))
    Y_f = np.insert(Y_f,0,Y[0])


    # transform back the signal to time domain
    iY = np.fft.ifft(Y_f)


    print('gc ' + str(len(gc)))
    print('len(y) ' + str(len(y)))
    print('len(iY) ' + str(len(iY)))
    print('len(Y_f) ' + str(len(Y_f)))
    print('len(Yp) ' + str(len(Yp)))
    print('len(Yp_f) ' + str(len(Yp_f)))
    print('len(t) ' + str(len(t)))
    

    # plot comparison of input and filtered by applying transfer function directly to freq. spectrum of signal
    fig4, ax4 = plt.subplots(2,1)
    ax4[0].plot(t,y, color='black') 
    ax4[0].plot(t,iY, color='red') 
    ax4[0].set_xlabel('Time')
    ax4[0].set_ylabel('Arb. Units.')
    ax4[0].set_title('Transfer function filtered response')

    ax4[1].plot(frq,np.abs(Yp), color='black') 
    iYp = iY[range(1,(n-1)/2+1)] # one sided
    #ax4[1].plot(frq,np.abs(iYp), color='red') # amplitude freq.
    ax4[1].plot(frq,np.abs(Yp_f), color='red') # amplitude freq.
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
    print('t_y_c_samples ' + str(t_y_c_samples))
    n_samples = len(y_c_samples)
    print(' n_samples ' + str( n_samples) )
    # do a FFT on sampled signal
    Y_c_samples = np.fft.fft(y_c_samples)
    Fs_samples = Fs*Fs_frac
    T_c_samples = len(y_c_samples)/Fs_samples # length in time of signal
    frq_samples = np.arange(len(y_c_samples))/T_c_samples # freq axis
    print('Fs_samples ' + str(Fs_samples))
    print('T_c_samples ' + str(T_c_samples))

    if n_samples%2==0:
        max_pos = len(Y_c_samples)/2    
        max_neg = len(Y_c_samples)/2
    else:
        max_pos = (len(Y_c_samples)-1)/2
        max_neg = (len(Y_c_samples)+1)/2

    Yp_c_samples = Y_c_samples[range(1,max_pos+1)]
    Ym_c_samples = Y_c_samples[range(max_neg,len(Y_c_samples))]
    # one-sided
    frq_samples = frq_samples[range(1,max_pos+1)]

    print(' frq_samples.shape ' + str( frq_samples.shape) )
    print(' frq_samples ' + str( frq_samples) )
    print(' Y_c_samples.shape ' + str( Y_c_samples.shape) )
    print(' Y_c_samples ' + str( Y_c_samples) )
    print(' Yp_c_samples.shape ' + str( Yp_c_samples.shape) )
    print(' Yp_c_samples ' + str( Yp_c_samples) )
    print(' Ym_c_samples.shape ' + str( Ym_c_samples.shape) )
    print(' Ym_c_samples ' + str( Ym_c_samples) )

    # deconvolute with filter transfer function
    #K_samples = 2*math.pi*frq_samples*tau     # helper constant, use same freq range as signal fft to apply
    if args.filter == 'gaus':
        gc_samples = tf(frq_samples, 1.0/(2*np.pi*tau))
    elif args.filter == 'rc':
        gc_samples = tf(frq_samples, tau)

    print(' gc_samples.shape ' + str( gc_samples.shape) )
    print(' gc_samples ' + str( gc_samples) )

    #gc_samples = gc_samples*np.sum(Yp_c_samples)/np.sum(gc_samples)

    print(' gc_samples ' + str( gc_samples) )

    #1/np.sqrt(1 + np.square(K_samples))     # direct transfer function for filter
    Yp_c_samples_deconv = Yp_c_samples/gc_samples # deconvolve
    Ym_c_samples_deconv = np.flipud(np.flipud(Ym_c_samples)/gc_samples)

    print(' Yp_c_samples_deconv ' + str( Yp_c_samples_deconv) )
    #print(' sum(Yp_c_samples_deconv) ' + str( np.sum(Yp_c_samples_deconv) ))

    #Ym_c_samples_deconv[1:len(Yp_c_samples_deconv)] = np.flipud(Yp_c_samples_deconv[1:len(Yp_c_samples_deconv)])
    #Yp_c_samples_deconv = Yp_c_samples_deconv*np.sum(Yp_c_samples)/np.sum((Yp_c_samples_deconv))
    #Ym_c_samples_deconv = Ym_c_samples_deconv*np.sum(Ym_c_samples)/np.sum((Ym_c_samples_deconv))    

    if n_samples%2==0:
        Y_c_samples_deconv = np.concatenate((Yp_c_samples_deconv,Ym_c_samples_deconv[1:len(Ym_c_samples_deconv)])) # max frequency at the same bin ?
    else:
        Y_c_samples_deconv = np.concatenate((Yp_c_samples_deconv,Ym_c_samples_deconv))
    Y_c_samples_deconv = np.insert(Y_c_samples_deconv,0,Y_c_samples[0])

    print(' Yp_c_samples_deconv ' + str( Yp_c_samples_deconv) )
    print(' Ym_c_samples_deconv ' + str( Ym_c_samples_deconv) )

    # inverse FFT to time domain    
    y_c_samples_deconv = np.fft.ifft(Y_c_samples_deconv)
    y_c_samples_deconv = y_c_samples_deconv*np.sum(y_c_samples)/np.sum((y_c_samples_deconv)) # normalize



    print(' len(Y_c_samples) ' + str( len(Y_c_samples)))
    print(' Fs ' + str( Fs ))
    print(' Fs_frac ' + str( Fs_frac) )
    print(' Fs_samples ' + str( Fs_samples) )
   # print(' frq_c ' + str( frq_samples) )
    #print(' frq_samples ' + str( frq_samples) )
    print(' y_c_samples.shape ' + str( y_c_samples.shape))
    print(' Y_c_samples.shape ' + str( Y_c_samples.shape))
    print(' Yp_c_samples.shape ' + str( Yp_c_samples.shape) )
    print(' Ym_c_samples.shape ' + str( Ym_c_samples.shape) )
    print(' Yp_c_samples_deconv.shape ' + str( Yp_c_samples_deconv.shape) )
    print(' Ym_c_samples_deconv.shape ' + str( Ym_c_samples_deconv.shape) )
    print(' Y_c_samples_deconv.shape ' + str( Y_c_samples_deconv.shape) )
    print(' y_c_samples_deconv.shape ' + str( y_c_samples_deconv.shape) )
    print(' y_c_samples_deconv ' + str( y_c_samples_deconv) )


    
    # plot comparison of input and filtered using convolution
    fig5, ax5 = plt.subplots(2,1)
    fig5.set_size_inches(12,6)
    ax5[0].plot(t,y, color='black', label='Input signal') 
    ax5[0].plot(t_y_c,y_c, color='red', label='Filtered') 
    ax5[0].plot(t_y_c_samples,y_c_samples, color='blue', linestyle=':', marker='o', label='Samples') 
    ax5[0].plot(t_y_c_samples,np.abs(y_c_samples_deconv), color='green', label='Deconvoluted')
    ax5[0].set_xlabel('Time')
    ax5[0].set_ylabel('Arb. Units.')
    ax5[0].set_title('ADC Response')
    legend = ax5[0].legend(loc='upper right', shadow=False)

    print('np.sum(y_c_samples) = ' + str(np.sum(y_c_samples)))
    print('np.sum(y_c_samples_deconv) = ' + str(np.sum(np.abs(y_c_samples_deconv))))

    
    #ax5[1].semilogy(frq_samples, np.abs(Y_c_samples), color='blue') 
    ax5[1].semilogy(frq_samples, np.abs(Yp_c_samples), color='blue',label='Yp_c_samples') 
    ax5[1].semilogy(frq_samples, np.abs(np.flipud(Ym_c_samples)), color='blue', label='Ym_c_samples') 
    ax5[1].semilogy(frq_samples, np.abs(Yp_c_samples_deconv), color='green', label='Yp_c_samples_deconv')  
    ax5[1].semilogy(frq_samples, np.abs(np.flipud(Ym_c_samples_deconv)), color='yellow',label='Ym_c_samples_deconv') 
   #ax5[1].semilogyy(frq_samples, np.abs(np.flipud(Ym_c_samples_deconv)), color='black') 
    #ax5[1].semilogy(np.abs(Y_c_samples_deconv), color='black') 
    ax5[1].semilogy(frq_samples, np.abs(gc_samples), color='cyan',label='gc_samples') 
    #ax5[1].semilogy(frq_samples[range(len(y_c_samples)/2)],gc_samples[range(len(y_c_samples)/2)], color='red') 
    #ax5[1].semilogy(frq_samples[range(len(y_c_samples)/2)],np.abs(Yp_c_samples_deconv), color='black') 
    ax5[1].set_xlabel('Frequency [Hz]')
    ax5[1].set_ylabel('Arb. Units.')
    legend = ax5[1].legend(loc='upper right', shadow=True)

    
    ax5[0].text(0.01, 0.9, os.path.basename(args.file[0]), transform=ax5[0].transAxes, color='green', fontsize=10)
    ax5[0].text(0.01, 0.8, r'$\tau$={0:.2f}$\mu$s'.format(tau*1e6), transform=ax5[0].transAxes, color='green', fontsize=10)
    ax5[0].text(0.01, 0.7, 'ADC {0:.2f}Msps'.format(args.freq), transform=ax5[0].transAxes, color='green', fontsize=10)
    

    fig5.savefig(saveFilename + '_filtered_signal_sampled.png',bbox_inches='tight')

    #ans = raw_input('dd')


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




    fig0, ax0 = plt.subplots(4,1)
    ax0[0].semilogy(frq_samples,1/abs(gc_samples)) # time signal
    ax0[1].plot(t_ir,y_ir) # impulse response
    yc_tmp = np.convolve(y,y_ir)
    ax0[2].plot(yc_tmp) # impulse response
    yc_tmp2 = np.convolve(y_ir,y)
    ax0[2].plot(yc_tmp2,color='red') # impulse response
    ax0[3].plot(yc_tmp/yc_tmp2,color='black') # impulse response
    print('y shape ' + str(y.shape))
    print('y_ir shape ' + str(y_ir.shape))
    print('yc_tmp shape ' + str(yc_tmp.shape))
    print('yc_tmp2 shape ' + str(yc_tmp2.shape))


    if not args.batch:
        plt.show()

    print('Done processing' + args.file[0])

