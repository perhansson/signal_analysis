
import sys
import re
import os
import pylab as plt
import math
import numpy as np
import argparse
#import matplotlib.pyplot as plt
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

    vals_max = []
    vals_axis = []
    for filename in args.file:
        print(filename)
        m = re.match('.*track_(\w)_d.*',filename)
        #print(m.group(1))
        d = m.group(1)
        #print(d)
        f = np.load(filename)
        #print(f.files)

        
        vals_max.append(np.max(f['yd']))
        vals_axis.append(d)
        #vals_min.append(np.min(f['yd']))

    x = range(len(args.file))
    ax1.plot(x,vals_max)
    plt.xticks(x,vals_axis,rotation='horizontal')
    ax1.set_xticklabels(vals_axis) 
    print(ax1.get_xticks().tolist())
    #ax1.set_xlabel('Type')
    ax1.set_ylabel('Max amplitude')
    fig1.tight_layout()

    if args.name != '':
        fig1.savefig(args.name)
    plt.show()
