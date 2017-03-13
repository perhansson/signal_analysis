
import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar



if __name__ == '__main__':
    print('Just GO')

    filename = sys.argv[1]
    tau = float(sys.argv[2])*1e-6


    Ts = 0.1e-6 # sampling interval
    Fs = 1.0/Ts # sampling rate
    

    # find histogram of time signal
    
    a = np.loadtxt(filename)
    #print(a)
    y = a[:,1]
    t = a[:,0]*1e-6


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

    #print(range((n-1)/2+1))
    #print(range((n+1)/2,n))    
    #print('n='+str(n))
    #print('Yp= ' + str(Yp))
    #print(Yp.shape)
    #print(Ym.shape)
    #ans = raw_input('continue?')


    # plot time signal and frequency spectrum of input signal
    fig1, ax1 = plt.subplots(2,1)
    ax1[0].plot(t,y) # time signal
    ax1[0].set_xlabel('Time (microsec)')
    ax1[0].set_ylabel('Current response (microA)')
    ax1[1].plot(frq,np.abs(Yp),'r') #color='red') # amplitude freq.
    ax1[1].plot(frq,np.flipud(np.abs(Ym)),'c') #color='red') # amplitude freq.
    ax1[1].set_xlabel('Frequency spectrum [Hz]')
    ax1[1].set_ylabel('Amplitude')

    fig1.savefig(os.path.splitext(os.path.basename(filename))[0] + '_input_signal.png',bbox_inches='tight')
    

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
    ax2[1].plot(frq_ir,gc_ir) 
    #ax2[1].semilogx(frq,gc) # freq
    #ax2[1].semilogx(frq,np.abs(Y_ir)) # freq
    ax2[1].plot([fc,fc],[0,1], color='k', linestyle='-', linewidth=2) 
    ax2[1].plot([fc/10,fc*10],[1/math.sqrt(2),1/math.sqrt(2)], color='k', linestyle='-', linewidth=2) 


    
    # convolve impulse respone with input signal
    print('y shape' + str(y.shape))
    print('y_ir shape' + str(y_ir.shape))    
    #y_c = np.convolve(y,y_ir)
    y_c = np.convolve(y_ir,y)
    y_c = y_c*np.sum(y)/np.sum(y_c)
    Y_C = np.fft.fft(y_c)
    Yp_C = Y_C[range(len(Y_C)/2)] # one sided

    print('y_ir' + str(y_ir))
    print('y' + str(y))
    print('y_c' + str(y_c))

    print('y_c shape' + str(y_c.shape))    
    print('Yp_C shape' + str(Yp_C.shape))    

    print('fc ' + str(fc))
    print('t_ir ' + str(t_ir))
    print('frq ' + str(frq))
    #frq = frq[range(n/2)]
    print('gc ' + str(gc))
    print(Yp.shape)
    print(frq.shape)
    print(gc.shape)

    
    # apply transfer function to frequency domain signal repr. 
    Yp_f = Yp*gc 
    Ym_f = np.flipud(np.flipud(Ym)*gc) 
    print(Yp_f.shape)
    Y_f = np.concatenate((Yp_f,Ym_f))
    Y_f = np.insert(Y_f,0,Y[0])
    print('Yp_f ' + str(Yp_f))
    print('Ym_f ' + str(Ym_f))
    print('Y_f ' + str(Y_f))


    # time domain tranformation
    iY = np.fft.ifft(Y_f)
    print('iY ' + str(iY))
    print(Yp.shape)
    print(Ym.shape)
    print(Yp_f.shape)
    print(Ym_f.shape)
    print(Y_f.shape)
    print(iY.shape)
    print(t.shape)




    # plot comparison of input and gc filtered
    fig3, ax3 = plt.subplots(2,1)
    #t_y_c = range(len(y_c))    
    print('len(y_c) = ' + str(len(y_c)) + ' len(t) = ' + str(len(t)) + ' n = ' + str(n))
    print('t = ' + str(t))
    print('Ts = ' + str(Ts))
    c = (np.arange(len(y_c)-n)*Ts + Ts + t[-1])
    print('c = ' + str(c))
    t_y_c = np.concatenate((t,c))
    T_c  = len(t_y_c)*Ts # length of signal
    frq_c = np.arange(len(t_y_c))/T_c


    print('Ts = ' + str(Ts) + ' Fs = ' + str(Fs))
    print('len(t_y_c) = ' + str(len(t_y_c)))
    print('t_y_c = ' + str(t_y_c))
    ax3[0].plot(t,y, color='black') 
    ax3[0].plot(t_y_c,y_c, color='red') 
    ax3[1].plot(frq,np.abs(Yp), color='black') 
    ax3[1].plot(frq_c[range(len(t_y_c)/2)],np.abs(Yp_C), color='red') 
                           
    #ans  = raw_input('c')


    # plot comparison of input and convolution filtered
    fig4, ax4 = plt.subplots(2,1)
    ax4[0].plot(t,y, color='black') 
    ax4[0].plot(t,iY, color='red') 
    ax4[1].plot(frq,np.abs(Yp), color='black') 
    iYp = iY[range(1,(n-1)/2+1)] # one sided
    ax4[1].plot(frq,np.abs(iYp), color='red') #color='red') # amplitude freq.





    print('sum y ' + str(np.sum(y)))
    print('sum y_ir ' + str(np.sum(y_ir)))
    print('sum y_c ' + str(np.sum(y_c)))

    




    plt.show()

    #ans = raw_input('continue?')
