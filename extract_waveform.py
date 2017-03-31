
import sys
import os
from ROOT import TFile

    
def process(name):
    # find histogram of time signal
    rootFilename = sys.argv[1]
    tFile = TFile(rootFilename)
    #h = tFile.Get('FieldResponse_Y')
    h = tFile.Get(name)
    dirname = os.path.dirname(rootFilename)
    plane = name.split('_')[1]
    asciiname = os.path.join(dirname, os.path.splitext(os.path.basename(rootFilename))[0] + '_' + plane + '.txt')
    print('Write: ' + asciiname)
    with open(asciiname,'w') as f:
        for i in range(1,h.GetNbinsX()):
            y = h.GetBinContent(i)
            x = h.GetBinCenter(i)
            s = '{0:.5f} {1:.14f}'.format(x,y) + '\n'
            f.write(s)


if __name__ == '__main__':
    print('Just GO')
    for name in ['FieldResponse_Y', 'FieldResponse_U','FieldResponse_V']:
        process(name)

