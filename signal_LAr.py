import ROOT
import sys
import math
import numpy as np
import matplotlib.pyplot as plt


# ionization energy in LAr
W_Lar = 23.6 #eV per e-h pair
# cell size
dx_cell = 0.5 #cm 
# MPV for energy loss in cell (Landau)
mpv_eloss = 0.8e6 #eV in 0.5cm 
#recombination for MIP
recomb = 0.7
# electron lifetime
lifetime = 3e-3 #seconds
# drift velocity 
v_drift = 0.15e6 #cm/second (assumes 500V/cm)
q_e = 1.6e-19 #C
e_fC = 1.0e-15/q_e # electrons per fC
rho_Lar = 1.3954 #g/cm3
N_A = 6.022e23
# noise
noise = 500.0*W_Lar #eV

# find MPV distribution for 
print('GO')

def mpv_charge(x, dedx, dx, recomb, t):
    # x: drift distance in cm
    # dedx: Eloss per cm in eV
    # dx: length of track
    # recomb: recombination factor
    # lt: lifetime
    td = x/v_drift
    print('td/t ' + str(td/t))
    print('1.0/exp ' + str(1.0/np.exp(td/t)))
    r = dedx*dx*recomb*1.0/np.exp(td/t)
    print(r)
    return r


# drift distance
x_drift = np.arange(360)
t_drift = (x_drift/v_drift) #second

#relative loss from lifetime
loss_x = 1.0/np.exp(t_drift/lifetime)

# MPV charge
q_mpv = mpv_charge(x_drift, mpv_eloss/dx_cell, dx_cell, recomb, lifetime) #eV

# plot
fig1,ax1 = plt.subplots()
fig1.patch.set_facecolor('white')
#axarr[1,0].set_bgcolor('grey
ax1.set_axis_bgcolor('white')
ax1.plot(x_drift,loss_x)
ax1.set_xlabel('Drift distance (cm)')
ax1.set_ylabel('Relative loss from lifetime ({0:.1f} ms)'.format(lifetime*1e3))
fig1.tight_layout()
fig1.savefig('q_mpv_loss_time_dune.png')

fig2, ax2 = plt.subplots()
fig2.patch.set_facecolor('white')
ax2.plot(x_drift,q_mpv/W_Lar)
ax2.set_ylabel('MPV Charge (e-)')
ax2.set_xlabel('Drift distance (cm)')
ax22 = ax2.twinx()
ax22.plot(x_drift, q_mpv/W_Lar/e_fC,color='red')
ax22.set_ylabel('MPV Charge (fC)',color='red')
ax22.tick_params('y', colors='r')
#fig2.tight_layout()
fig2.savefig('q_mpv_dune.png')



SN_req = 9.0
N_req = (q_mpv/W_Lar)/SN_req

fig3, ax3 = plt.subplots()
fig3.patch.set_facecolor('white')
ax3.plot(x_drift,N_req/1.5, 'k', label='0deg', color='black')
ax3.plot(x_drift,N_req/1.5/math.sin(math.pi/180*45), 'k:', label='45deg (S/sin(45))',color='black')
ax3.plot(x_drift,N_req/4.8/math.sin(math.pi/180*45), 'k--', label='45deg (S/sin(45)/4.8)',color='black')
ax3.set_ylabel('Max noise (e-)')
ax3.set_xlabel('Drift distance (cm)')
ax3.legend(loc='upper right')
ax32 = ax3.twinx()
ax32.plot(x_drift,N_req/1.5/e_fC, 'k', label='0deg', color='red')
ax32.plot(x_drift,N_req/1.5/math.sin(math.pi/180*45)/e_fC, 'k:', label='45deg (S/sin(45))',color='red')
ax32.plot(x_drift,N_req/4.8/math.sin(math.pi/180*45)/e_fC, 'k--', label='45deg (S/sin(45)/4.8)',color='red')
ax32.set_ylabel('Max noise (fC)',color='red')
ax32.tick_params('y', colors='r')
#fig3.tight_layout()
fig3.savefig('noise_mpv_dune.png')




plt.show()


#dedx = mpv_eloss/dx_cell
#q_mip = mpv_charge(x, dedx, recomb, lifetime)


#def xsi2(z,Z,A,beta,rho):
#    KoverA = 0.31*A
#    return KoverA/2*Z/A*1/(beta*beta)
#def xsi(z,Z,A,beta,rho):
#    return 0.1535*z*z*Z/(A*beta*beta)*rho


#xsi_L_si = xsi2(1.0, 14 , 28, 1.0, 2.3) 
#w_L = 4.0*dx_cell*xsi(1.0, 18, 40, 1.0, rho_Lar) 

#print('xsi_L ' + str(xsi_L_si))
#print('w_L ' + str(w_L))
#print('w_L_si ' + str(w_L_si/1000.0) + ' keV')

l0 = ROOT.TF1("l0","TMath::Landau(x,0)",-1,5) #,-mpv_eloss/W_Lar,mpv_eloss/W_Lar*5)
#l0.SetParameters(mpv_eloss/W_Lar)


c0 = ROOT.TCanvas('c0','c0',10,10,700,500)
l0.Draw()

ans = raw_input('continue?')


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


