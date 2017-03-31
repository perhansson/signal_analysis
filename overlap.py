import ROOT
import sys
import math



mpv = float(sys.argv[1])
noise = float(sys.argv[2])


print('GO')
print('MPV ' + str(mpv) + ' noise ' + str(noise))


l1 = ROOT.TF1("l1","TMath::Landau(x,[0],[1],0)",0,100000);
l1.SetParameters(mpv,noise);

l2 = ROOT.TF1("l2","TMath::Landau(x,[0],[1],0)",0,100000);
l2.SetParameters(mpv*2,noise);


c1 = ROOT.TCanvas('c1','c1',10,10,700,500)
l1.Draw()
l2.Draw("L, same") 

tot2 = l2.Integral(0,100000)
print('tot2 = ' + str(tot2))

s2_thr = 0.9
a2_low =100000
r2 = 0.0

while r2 < s2_thr:
    a2 = l2.Integral(a2_low,100000)
    r2 = a2/tot2
    print('a2_low ' + str(a2_low) + ' a2 ' + str(a2) + ' tot2 ' + str(tot2) + ' r2 ' + str(r2))
    a2_low -= 100.0

print('a2_low ' + str(a2_low) + ' -> r2 ' + str(r2))


tot1 = l1.Integral(0,100000)
r1 = 1 - l1.Integral(0,a2_low)/tot1

print('a2_low ' + str(a2_low) + ' -> r1 ' + str(r1))






print('\n\n\n')


def ms_angle(E,p,z,x,X0):

    beta = p/E
    c = 1.0
    return 13.6/(beta*p*c)*z*math.sqrt(x/X0)*(1.0 + 0.038*math.log(x/X0))    


E = 1000.0
p = 1000.0
E = math.sqrt(p*p + 105.0*105.0)
Z = 1.0
X0 = 14.00 #cm
y_rms = ROOT.TGraph()
i = 0
for ix in range(1,101):
    x = ix * 0.1
    th = ms_angle(E,p,Z,x,X0)
    y = 1.0/math.sqrt(3.0)*th*x*X0 #cm
    print('x ' + str(x) + '/' + str(x*X0) + 'cm th ' + str(th) + ' y ' + str(y))
    y_rms.SetPoint(i,x*X0,y)
    i+=1

c2 = ROOT.TCanvas('c2','c2',50,50,750,550)
y_rms.Draw('ALP')
y_rms.SetTitle('MS displacement p=1Gev/c muon;Distance travelled (cm);Displacement(cm)')







ans = raw_input('cont.?')


