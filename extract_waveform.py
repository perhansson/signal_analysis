
import sys
import os
from ROOT import TFile

if __name__ == '__main__':
    print('Just GO')
    # find histogram of time signal
    rootFilename = sys.argv[1]
    tFile = TFile(rootFilename)
    h = tFile.Get('FieldResponse_Y')
    dirname = os.path.dirname(rootFilename)
    asciiname = os.path.join(dirname, os.path.splitext(os.path.basename(rootFilename))[0] + '.txt')
    print('Write: ' + asciiname)
    with open(asciiname,'w') as f:
        for i in range(1,h.GetNbinsX()):
            y = h.GetBinContent(i)
            x = h.GetBinCenter(i)
            s = '{0:.5f} {1:.14f}'.format(x,y) + '\n'
            f.write(s)


