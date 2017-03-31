
import sys
import re
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
    parser.add_argument('name', help='output name.')
    parser.add_argument('file', nargs='*', help='Data files.')
    args = parser.parse_args()
    print(args)
    return args



if __name__ == '__main__':
    print('Just GO')

    args = get_args()

    fig1,ax1 = plt.subplots()
    fig1.patch.set_facecolor('white')
    #axarr[1,0].set_bgcolor('grey
    ax1.set_axis_bgcolor('white')

    

    for filename in args.file:
        print(filename)
        m = re.match('.*_d(\d.*)cm.*',filename)
        #print(m.group(1))
        d = int(float(m.group(1)))
        #print(d)
        f = np.load(filename)
        #print(f.files)
        ax1.plot(f['t_ext2']*1e6,f['yd'])

        print('max ' + str(np.max(f['yd'])))
        print('min ' + str(np.min(f['yd'])))

        sum_p = np.sum(f['yd'][f['yd']>0])
        sum_m = np.sum(f['yd'][f['yd']<=0])
        sum_all = np.sum(f['yd'])
        print('sum_p ' + str(sum_p))
        print('sum_m ' + str(sum_m))
        print('sum ' + str(sum_all))
    
        
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Signal')
    fig1.tight_layout()

    if args.name != '':
        fig1.savefig(args.name)
    plt.show()

            
        
