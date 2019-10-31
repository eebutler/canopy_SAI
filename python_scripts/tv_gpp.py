#cd /Users/timaeus/Desktop/Projects/Papers_InProgress/SAI_canopy/output

import numpy as np
import matplotlib.pyplot as plt

tveg = np.loadtxt('avg_tveg.txt')
gpp = np.loadtxt('avg_gpp.txt')

# transform to vectors
tveg = np.ravel(tveg)
gpp = np.ravel(gpp)

# update missing values with NaN
tveg[tveg==-999] = np.nan
gpp[gpp==-999] = np.nan

# divide into "modest" and "substantial" T change
#idx = np.where((tveg<0.5)&(tveg>-0.5))
#idx2 = np.where((tveg>0.5)|(tveg<-0.5))


#plt.plot(tveg[idx],gpp[idx],'b.',label='|T$_{diff}|<0.5^{\circ}$C')
#plt.plot(tveg[idx2],gpp[idx2],'r.',label='|T$_{diff}|>0.5^{\circ}$C')
plt.plot(tveg,gpp,'k.')
plt.xlabel('Annual Average Vegetation Temperature Difference [$^{\circ}$C]')
plt.ylabel('Annual Average GPP Difference [gC m$^{-2}$ yr$^{-1}$]')
#plt.legend()
plt.savefig('tv_gpp.pdf')

# correlations
np.corrcoef(tveg[~np.isnan(tveg)],gpp[~np.isnan(gpp)])
#np.corrcoef(tveg[idx],gpp[idx])
#np.corrcoef(tveg[idx2],gpp[idx2])
