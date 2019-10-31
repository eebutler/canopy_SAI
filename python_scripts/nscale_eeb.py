# cd Desktop/Projects/Papers_InProgress/SAI_canopy/

import numpy as np
import matplotlib.pyplot as plt

# nitrogen extinction, contested values see Bonan et al 2011 and 2012
kn = 0.3 
# kb (empirical optical depth), 1.5 is approximately 70 deg. sun angle
# from Ross (1975), Goudriaan (1977), and Dai, et al. (2004)
chi = 0 # spherical leaf distribution, needleleaf trees = 0.01
zen = 70.0
cosz = np.cos(zen*np.pi/180)
phi1 = 0.5 - 0.633*chi - 0.33*chi**2
phi2 = 0.877 * (1 - 2*phi1)
gdir = phi1 + phi2*cosz
kb = gdir/cosz
# kb = 1.5

# singularity at lai = 0 so we omit that point, fixed sai
lai = np.arange(0.01,6.01,0.01)

def Kb(zen,chi=0):
    cosz = np.cos(zen*np.pi/180)
    phi1 = 0.5 - 0.633*chi - 0.33*chi**2
    phi2 = 0.877 * (1 - 2*phi1)
    gdir = phi1 + phi2*cosz
    kb = gdir/cosz
    return kb

def fsun(lai,sai,kb):
    fsu = (1 - np.exp(-kb*(lai+sai))) / (kb*(lai+sai))
    return fsu

def na_osun(lai,kb,kn):
    na = (1 - np.exp(-(kn+kb)*lai)) / (kn+kb)
    return na

#chi = 0
Fsun1 = na_osun(lai,Kb(70,chi),0.3)/(fsun(lai,0.5,Kb(70,chi))*lai)
Fsun2 = na_osun(lai,Kb(0,chi),0.3)/(fsun(lai,5,Kb(0,chi))*lai)
Fsun3 = na_osun(lai,Kb(70,chi),0.3)/(fsun(lai,5,Kb(70,chi))*lai)
Fsun4 = na_osun(lai,Kb(0,chi),0.3)/(fsun(lai,0.5,Kb(0,chi))*lai)

# Figure 1
# uncomment plt.plot below for the other calculations with equivalent SAI and zenith
plt.plot(lai,[1]*len(lai),'k')
plt.plot(lai,Fsun1,'g--',label='zenith = 70$^{\circ}$, SAI = 0.5')
plt.plot(lai,Fsun2,':',color='brown',label='zenith = 0$^{\circ}$  , SAI = 5.0')
#plt.plot(lai,Fsun3,'--',color='brown',label='zenith = 70$^{\circ}$, SAI = 5.0')
#plt.plot(lai,Fsun4,'g:',label='zenith = 0$^{\circ}$  , SAI = 0.5')
plt.xlabel('LAI')
plt.ylabel('Sunlit $N_{a}$')
plt.legend()
plt.savefig('na_old.pdf')
plt.close('all')

# Fixed SAI 5
sai = 0.5 

## old Fsun
Fsun = (1 - np.exp(-kb*(lai+sai))) / (kb*(lai+sai))
# nscaler in sun and shade
nscaler_osun = (1 - np.exp(-(kn+kb)*lai)) / (kn+kb)
nscaler_osha = (1 - np.exp(-kn*lai)) / kn - nscaler_osun
# lai scaled nscaler
Nscaler_osun = nscaler_osun/(Fsun*lai)
Nscaler_osha = nscaler_osha/((1-Fsun)*lai)

## new Fsun
Fsun_new = (1-np.exp(-kb*lai))*(1-np.exp(-kb*sai))/(kb*kb*sai*lai)

Nscaler_sun = ( ((1-np.exp(-kb*sai))/kb) * (1-np.exp(-(kn+kb)*lai))/(kn+kb) )/(Fsun_new*lai*sai) 
Nscaler_sha = ( sai*( 1-np.exp(-kn*lai) )/kn - \
              ( ((1-np.exp(-kb*sai))/kb) * (1-np.exp(-(kn+kb)*lai))/(kn+kb)) )/ \
              ( (1 - Fsun_new) * lai * sai ) 

# whole canopy nscaler (no optical extinction)
nscaler_whole = (1 - np.exp(-kn*lai)) / kn
Nscaler_whole = nscaler_whole/lai


# Figure 2
plt.plot(lai,Fsun,'g--',label='additive model')
plt.plot(lai,Fsun_new,'g',label='multiplicative model')
plt.xlabel('LAI')
plt.ylabel('Sunlit Fraction of Canopy')
plt.legend()
plt.xlim([0,6])
plt.ylim([0,1])
plt.savefig('fsun.pdf')
plt.close('all')

# Figure 3
plt.plot(lai,Nscaler_osun,'g--',label='old sun')
plt.plot(lai,Nscaler_osha,'b--',label='old shade')
plt.plot(lai,Nscaler_sun,'g',label='new sun')
plt.plot(lai,Nscaler_sha,'b',label='new shade')
plt.plot(lai,Nscaler_whole,'r',label='whole canopy')
plt.xlabel('LAI')
plt.ylabel('N$_{a}$')
plt.legend()
plt.xlim([0,6])
plt.ylim([0,1.5])
plt.savefig('na_compare.pdf')
plt.close('all')


#### a little bit more to play with not currently displayed in manuscript
def fsun_new(lai,sai,kb):
    fsu_new = (1-np.exp(-kb*lai))*(1-np.exp(-kb*sai))/(kb*kb*sai*lai)
    return fsu_new

def na_sun(lai,sai,kb,kn):
    na = ( ((1-np.exp(-kb*sai))/kb) * (1-np.exp(-(kn+kb)*lai))/(kn+kb) )
    return na

chi = 0
Fsun1 = na_osun(lai,Kb(70,chi),0.3)/(fsun(lai,0.5,Kb(70,chi))*lai)
Fsun2 = na_osun(lai,Kb(0,chi),0.3)/(fsun(lai,5,Kb(0,chi))*lai)
Fsun3 = na_sun(lai,0.5,Kb(70,chi),0.3)/(fsun_new(lai,0.5,Kb(70,chi))*lai*0.5)
Fsun4 = na_sun(lai,5.0,Kb(0,chi),0.3)/(fsun_new(lai,5,Kb(0,chi))*lai*5)

plt.plot(lai,Fsun1,'g--',label='Additive, zenith = 70$^{\circ}$, SAI = 0.5')
plt.plot(lai,Fsun2,'g:',label='Additive, zenith = 0$^{\circ}$  , SAI = 5.0')
plt.plot(lai,Fsun3,'b--',label='Multiplicative, zenith = 70$^{\circ}$, SAI = 0.5')
plt.plot(lai,Fsun4,'b:',label='Multiplicative, zenith = 0$^{\circ}$  , SAI = 5.0')
plt.xlabel('LAI')
plt.ylabel('Sunlit N$_{a}$')
plt.legend()
#plt.show()
