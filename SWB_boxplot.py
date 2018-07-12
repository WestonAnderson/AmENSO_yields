# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 15:00:24 2016
@author: westonanderson

Edited to average the soil water content over a specific region for
either flowering season or the growing season (defined later in the script)
"""
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes
                
seas = 'flwr'
years = [1950,2011]
numYrs = years[1]-years[0]
ENyrs = np.array([1953, 1957, 1963, 1965, 1968, 1972, 1976, 1979, 1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009]) #EN
LNyrs = np.array([1954, 1964, 1970, 1973, 1983, 1988, 1995, 1998, 2007,2010]) #LN

ENind = np.repeat(np.array(ENyrs - years[0])*12,12) + np.tile(range(6,18),np.size(ENyrs))
ENyrInd = np.zeros([numYrs*12],dtype='int')
EN0ind =np.copy(ENyrInd); EN0ind[ENind] =1
EN_1ind =np.copy(ENyrInd); EN_1ind[ENind-12] =1
EN1ind =np.copy(ENyrInd); EN1ind[ENind[ENind<714]+12] =1
EN2ind =np.copy(ENyrInd); EN2ind[ENind[ENind<702]+24] =1
#NOTE: the extra '+12' (+24) is because an EN calendar year ends in dec, but soybean/maize season starts Jan of next calendar year


LNind = np.repeat(np.array(LNyrs - years[0])*12,12) + np.tile(range(6,18),np.size(LNyrs))
LNyrInd = np.zeros([numYrs*12],dtype='int')
LN0ind =np.copy(LNyrInd); LN0ind[LNind[LNind<732]] =1
LN_1ind =np.copy(LNyrInd); LN_1ind[LNind-12] =1
LN1ind =np.copy(LNyrInd); LN1ind[LNind[LNind<714]+12] =1
LN2ind =np.copy(LNyrInd); LN2ind[LNind[LNind<702]+24] =1

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#
# Calculate the soil water balance
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/project/AmENSO_yields/soilWaterBalance.py')
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/nan_helper.py')
#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#Flowering seasons for each Kc curve by crop. Each season is two months long
yrOffset = np.array(range((numYrs/3)+1))*36
if seas is 'grow':
    soy1 = [-2,-1,0,1,2,3,  10,11,12,13,14]
    soy2 = [-2,-1,0,1,2,  22,23,24,25,26,27]
    soy3 = [12,13,14,15,  22,23,24,25,26]
    maize1 = [22,23,24,25,26,27]
    maize2 = [10,11,12,13,14,15]
    maize3 = [-2,-1,0,1,2,3]
    safS = [8,9,10,11,12,  20,21,22,23,24,  32,33,34,35,36]
    safM = [1,2,3,4,  13,14,15,16,  25,26,27,28]
elif seas is 'flwr':
    soy1 = [1,2,  12,13]
    soy2 = [0,1,  25,26]
    soy3 = [13,14, 24,25]
    maize1 = [24,25]
    maize2 = [12,13]
    maize3 = [0,1]
    safS = [10,11,  22,23,  34,35]
    safM = [2,3,  14,15,  26,27]

Kc_S1 = np.tile(soy1,(numYrs/3)+1)+np.repeat(yrOffset,np.size(soy1))
#Kc_S1 = Kc_S1[6:-5]
Kc_S1=Kc_S1[Kc_S1<numYrs*12];Kc_S1=Kc_S1[Kc_S1>=0]
Kc_S2 = np.tile(soy2,(numYrs/3)+1)+np.repeat(yrOffset,np.size(soy2));
#Kc_S2 = Kc_S2[5:-6]
Kc_S2=Kc_S2[Kc_S2<numYrs*12];Kc_S2=Kc_S2[Kc_S2>=0]
Kc_S3 = np.tile(soy3,(numYrs/3)+1)+np.repeat(yrOffset,np.size(soy3));
#Kc_S3=Kc_S3[:-9]
Kc_S3=Kc_S3[Kc_S3<numYrs*12];Kc_S3=Kc_S3[Kc_S3>=0]

Kc_M1 = np.tile(maize1,(numYrs/3)+1)+np.repeat(yrOffset+12,np.size(maize1));
#Kc_M1 =Kc_M1[:-6] 
Kc_M1=Kc_M1[Kc_M1<numYrs*12];Kc_M1=Kc_M1[Kc_M1>=0]
Kc_M2 = np.tile(maize2,(numYrs/3)+1)+np.repeat(yrOffset+12,np.size(maize2));
#Kc_M2 = Kc_M2[:-6]
Kc_M2=Kc_M2[Kc_M2<numYrs*12];Kc_M2=Kc_M2[Kc_M2>=0]
Kc_M3 = np.tile(maize3,(numYrs/3)+1)+np.repeat(yrOffset+12,np.size(maize3));
#Kc_M3 = Kc_M3[:-6]
Kc_M3=Kc_M3[Kc_M3<numYrs*12];Kc_M3=Kc_M3[Kc_M3>=0]

Kc_saf_S = np.tile(safS,(numYrs/3)+1)+np.repeat(yrOffset,np.size(safS));
#Kc_saf_S=Kc_saf_S[:-15]
Kc_saf_S = Kc_saf_S[Kc_saf_S<numYrs*12];Kc_saf_S=Kc_saf_S[Kc_saf_S>=0]
Kc_saf_M = np.tile(safM,(numYrs/3)+1)+np.repeat(yrOffset+12,np.size(safM));
#Kc_saf_M=Kc_saf_M[:-12]
Kc_saf_M = Kc_saf_M[Kc_saf_M<numYrs*12];Kc_saf_M=Kc_saf_M[Kc_saf_M>=0]

s1_ind = np.zeros([numYrs*12],dtype='int'); s1_ind[Kc_S1]=1
s2_ind = np.zeros([numYrs*12],dtype='int'); s2_ind[Kc_S2]=1
s3_ind = np.zeros([numYrs*12],dtype='int'); s3_ind[Kc_S3]=1

m1_ind = np.zeros([numYrs*12],dtype='int'); m1_ind[Kc_M1]=1
m2_ind = np.zeros([numYrs*12],dtype='int'); m2_ind[Kc_M2]=1
m3_ind = np.zeros([numYrs*12],dtype='int'); m3_ind[Kc_M3]=1

s_saf_ind = np.zeros([numYrs*12],dtype='int'); s_saf_ind[Kc_saf_S]=1
m_saf_ind = np.zeros([numYrs*12],dtype='int'); m_saf_ind[Kc_saf_M]=1

#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~
######## Set the region of analysis ###############


lonMin = -57; lonMax = -47
latMin = -15; latMax = -25
            #-10
    
latMax = np.min(np.where(lats<latMax))
latMin = np.max(np.where(lats>latMin))
lonMin = np.max(np.where(lons<lonMin))
lonMax = np.min(np.where(lons>lonMax))


#general water stress timeseries.
ws1TS = np.nanmean(np.nanmean(w_s1[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)
ws2TS = np.nanmean(np.nanmean(w_s2[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)
ws3TS = np.nanmean(np.nanmean(w_s3[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)
wcerTS = np.nanmean(np.nanmean(w_cer[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)

#NOTE: the '+1' is because an EN calendar year ends in dec, 
#      but soybean/maize season starts jan or feb of next calendar year
s1_EN_1 = np.copy(ws1TS);s1_EN_1[~((s1_ind==1)&(EN0ind==1))]=np.nan;s1_EN_1=np.reshape(s1_EN_1[6:-6],[12,numYrs-1])
s1_LN_1 = np.copy(ws1TS);s1_LN_1[~((s1_ind==1)&(LN0ind==1))]=np.nan;s1_LN_1=np.reshape(s1_LN_1[6:-6],[12,numYrs-1])
s1_EN0 = np.copy(ws1TS);s1_EN0[~((s1_ind==1)&(EN1ind==1))]=np.nan;s1_EN0=np.reshape(s1_EN0[6:-6],[12,numYrs-1])
s1_LN0 = np.copy(ws1TS);s1_LN0[~((s1_ind==1)&(LN1ind==1))]=np.nan;s1_LN0=np.reshape(s1_LN0[6:-6],[12,numYrs-1])
s1_EN1 = np.copy(ws1TS);s1_EN1[~((s1_ind==1)&(EN2ind==1))]=np.nan;s1_EN1=np.reshape(s1_EN1[6:-6],[12,numYrs-1])
s1_LN1 = np.copy(ws1TS);s1_LN1[~((s1_ind==1)&(LN2ind==1))]=np.nan;s1_LN1=np.reshape(s1_LN1[6:-6],[12,numYrs-1])
s1_EN_1 = np.nanmean(s1_EN_1,0);s1_LN_1 = np.nanmean(s1_LN_1,0)
s1_EN0 = np.nanmean(s1_EN0,0);s1_LN0 = np.nanmean(s1_LN0,0)
s1_EN1 = np.nanmean(s1_EN1,0);s1_LN1 = np.nanmean(s1_LN1,0)

s2_EN_1 = np.copy(ws2TS);s2_EN_1[~((s2_ind==1)&(EN0ind==1))]=np.nan;s2_EN_1=np.reshape(s2_EN_1[6:-6],[12,numYrs-1])
s2_LN_1 = np.copy(ws2TS);s2_LN_1[~((s2_ind==1)&(LN0ind==1))]=np.nan;s2_LN_1=np.reshape(s2_LN_1[6:-6],[12,numYrs-1])
s2_EN0 = np.copy(ws2TS);s2_EN0[~((s2_ind==1)&(EN1ind==1))]=np.nan;s2_EN0=np.reshape(s2_EN0[6:-6],[12,numYrs-1])
s2_LN0 = np.copy(ws2TS);s2_LN0[~((s2_ind==1)&(LN1ind==1))]=np.nan;s2_LN0=np.reshape(s2_LN0[6:-6],[12,numYrs-1])
s2_EN1 = np.copy(ws2TS);s2_EN1[~((s2_ind==1)&(EN2ind==1))]=np.nan;s2_EN1=np.reshape(s2_EN1[6:-6],[12,numYrs-1])
s2_LN1 = np.copy(ws2TS);s2_LN1[~((s2_ind==1)&(LN2ind==1))]=np.nan;s2_LN1=np.reshape(s2_LN1[6:-6],[12,numYrs-1])
s2_EN_1 = np.nanmean(s2_EN_1,0);s2_LN_1 = np.nanmean(s2_LN_1,0)
s2_EN0 = np.nanmean(s2_EN0,0);s2_LN0 = np.nanmean(s2_LN0,0)
s2_EN1 = np.nanmean(s2_EN1,0);s2_LN1 = np.nanmean(s2_LN1,0)

s3_EN_1 = np.copy(ws3TS);s3_EN_1[~((s3_ind==1)&(EN0ind==1))]=np.nan;s3_EN_1=np.reshape(s3_EN_1[6:-6],[12,numYrs-1])
s3_LN_1 = np.copy(ws3TS);s3_LN_1[~((s3_ind==1)&(LN0ind==1))]=np.nan;s3_LN_1=np.reshape(s3_LN_1[6:-6],[12,numYrs-1])
s3_EN0 = np.copy(ws3TS);s3_EN0[~((s3_ind==1)&(EN1ind==1))]=np.nan;s3_EN0=np.reshape(s3_EN0[6:-6],[12,numYrs-1])
s3_LN0 = np.copy(ws3TS);s3_LN0[~((s3_ind==1)&(LN1ind==1))]=np.nan;s3_LN0=np.reshape(s3_LN0[6:-6],[12,numYrs-1])
s3_EN1 = np.copy(ws3TS);s3_EN1[~((s3_ind==1)&(EN2ind==1))]=np.nan;s3_EN1=np.reshape(s3_EN1[6:-6],[12,numYrs-1])
s3_LN1 = np.copy(ws3TS);s3_LN1[~((s3_ind==1)&(LN2ind==1))]=np.nan;s3_LN1=np.reshape(s3_LN1[6:-6],[12,numYrs-1])
s3_EN_1 = np.nanmean(s3_EN_1,0);s3_LN_1 = np.nanmean(s3_LN_1,0)
s3_EN0 = np.nanmean(s3_EN0,0);s3_LN0 = np.nanmean(s3_LN0,0)
s3_EN1 = np.nanmean(s3_EN1,0);s3_LN1 = np.nanmean(s3_LN1,0)

m1_EN_1 = np.copy(ws1TS);m1_EN_1[~((m1_ind==1)&(EN0ind==1))]=np.nan;m1_EN_1=np.reshape(m1_EN_1[6:-6],[12,numYrs-1])
m1_LN_1 = np.copy(ws1TS);m1_LN_1[~((m1_ind==1)&(LN0ind==1))]=np.nan;m1_LN_1=np.reshape(m1_LN_1[6:-6],[12,numYrs-1])
m1_EN0 = np.copy(ws1TS);m1_EN0[~((m1_ind==1)&(EN1ind==1))]=np.nan;m1_EN0=np.reshape(m1_EN0[6:-6],[12,numYrs-1])
m1_LN0 = np.copy(ws1TS);m1_LN0[~((m1_ind==1)&(LN1ind==1))]=np.nan;m1_LN0=np.reshape(m1_LN0[6:-6],[12,numYrs-1])
m1_EN1 = np.copy(ws1TS);m1_EN1[~((m1_ind==1)&(EN2ind==1))]=np.nan;m1_EN1=np.reshape(m1_EN1[6:-6],[12,numYrs-1])
m1_LN1 = np.copy(ws1TS);m1_LN1[~((m1_ind==1)&(LN2ind==1))]=np.nan;m1_LN1=np.reshape(m1_LN1[6:-6],[12,numYrs-1])
m1_EN_1 = np.nanmean(m1_EN_1,0);m1_LN_1 = np.nanmean(m1_LN_1,0)
m1_EN0 = np.nanmean(m1_EN0,0);m1_LN0 = np.nanmean(m1_LN0,0)
m1_EN1 = np.nanmean(m1_EN1,0);m1_LN1 = np.nanmean(m1_LN1,0)

m2_EN_1 = np.copy(ws2TS);m2_EN_1[~((m2_ind==1)&(EN0ind==1))]=np.nan;m2_EN_1=np.reshape(m2_EN_1[6:-6],[12,numYrs-1])
m2_LN_1 = np.copy(ws2TS);m2_LN_1[~((m2_ind==1)&(LN0ind==1))]=np.nan;m2_LN_1=np.reshape(m2_LN_1[6:-6],[12,numYrs-1])
m2_EN0 = np.copy(ws2TS);m2_EN0[~((m2_ind==1)&(EN1ind==1))]=np.nan;m2_EN0=np.reshape(m2_EN0[6:-6],[12,numYrs-1])
m2_LN0 = np.copy(ws2TS);m2_LN0[~((m2_ind==1)&(LN1ind==1))]=np.nan;m2_LN0=np.reshape(m2_LN0[6:-6],[12,numYrs-1])
m2_EN1 = np.copy(ws2TS);m2_EN1[~((m2_ind==1)&(EN2ind==1))]=np.nan;m2_EN1=np.reshape(m2_EN1[6:-6],[12,numYrs-1])
m2_LN1 = np.copy(ws2TS);m2_LN1[~((m2_ind==1)&(LN2ind==1))]=np.nan;m2_LN1=np.reshape(m2_LN1[6:-6],[12,numYrs-1])
m2_EN_1 = np.nanmean(m2_EN_1,0);m2_LN_1 = np.nanmean(m2_LN_1,0)
m2_EN0 = np.nanmean(m2_EN0,0);m2_LN0 = np.nanmean(m2_LN0,0)
m2_EN1 = np.nanmean(m2_EN1,0);m2_LN1 = np.nanmean(m2_LN1,0)

m3_EN_1 = np.copy(ws3TS);m3_EN_1[~((m3_ind==1)&(EN0ind==1))]=np.nan;m3_EN_1=np.reshape(m3_EN_1[6:-6],[12,numYrs-1])
m3_LN_1 = np.copy(ws3TS);m3_LN_1[~((m3_ind==1)&(LN0ind==1))]=np.nan;m3_LN_1=np.reshape(m3_LN_1[6:-6],[12,numYrs-1])
m3_EN0 = np.copy(ws3TS);m3_EN0[~((m3_ind==1)&(EN1ind==1))]=np.nan;m3_EN0=np.reshape(m3_EN0[6:-6],[12,numYrs-1])
m3_LN0 = np.copy(ws3TS);m3_LN0[~((m3_ind==1)&(LN1ind==1))]=np.nan;m3_LN0=np.reshape(m3_LN0[6:-6],[12,numYrs-1])
m3_EN1 = np.copy(ws3TS);m3_EN1[~((m3_ind==1)&(EN2ind==1))]=np.nan;m3_EN1=np.reshape(m3_EN1[6:-6],[12,numYrs-1])
m3_LN1 = np.copy(ws3TS);m3_LN1[~((m3_ind==1)&(LN2ind==1))]=np.nan;m3_LN1=np.reshape(m3_LN1[6:-6],[12,numYrs-1])
m3_EN_1 = np.nanmean(m3_EN_1,0);m3_LN_1 = np.nanmean(m3_LN_1,0)
m3_EN0 = np.nanmean(m3_EN0,0);m3_LN0 = np.nanmean(m3_LN0,0)
m3_EN1 = np.nanmean(m3_EN1,0);m3_LN1 = np.nanmean(m3_LN1,0)

m_saf_EN_1= np.copy(wcerTS); m_saf_EN_1[~((m_saf_ind==1)&(EN0ind==1))]=np.nan;m_saf_EN_1=np.reshape(m_saf_EN_1[6:-6],[12,numYrs-1])
m_saf_LN_1= np.copy(wcerTS); m_saf_LN_1[~((m_saf_ind==1)&(LN0ind==1))]=np.nan;m_saf_LN_1=np.reshape(m_saf_LN_1[6:-6],[12,numYrs-1])
m_saf_EN0= np.copy(wcerTS); m_saf_EN0[~((m_saf_ind==1)&(EN1ind==1))]=np.nan;m_saf_EN0=np.reshape(m_saf_EN0[6:-6],[12,numYrs-1])
m_saf_LN0= np.copy(wcerTS); m_saf_LN0[~((m_saf_ind==1)&(LN1ind==1))]=np.nan;m_saf_LN0=np.reshape(m_saf_LN0[6:-6],[12,numYrs-1])
m_saf_EN1= np.copy(wcerTS); m_saf_EN1[~((m_saf_ind==1)&(EN2ind==1))]=np.nan;m_saf_EN1=np.reshape(m_saf_EN1[6:-6],[12,numYrs-1])
m_saf_LN1= np.copy(wcerTS); m_saf_LN1[~((m_saf_ind==1)&(LN2ind==1))]=np.nan;m_saf_LN1=np.reshape(m_saf_LN1[6:-6],[12,numYrs-1])

#Safrinha soy is the only anomaly that happens in the same 'calendar year' as the enso event
s_saf_EN_1= np.copy(wcerTS); s_saf_EN_1[~((s_saf_ind==1)&(EN_1ind==1))]=np.nan;s_saf_EN_1=np.reshape(s_saf_EN_1[6:-6],[12,numYrs-1])
s_saf_LN_1= np.copy(wcerTS); s_saf_LN_1[~((s_saf_ind==1)&(LN_1ind==1))]=np.nan;s_saf_LN_1=np.reshape(s_saf_LN_1[6:-6],[12,numYrs-1])
s_saf_EN0= np.copy(wcerTS); s_saf_EN0[~((s_saf_ind==1)&(EN0ind==1))]=np.nan;s_saf_EN0=np.reshape(s_saf_EN0[6:-6],[12,numYrs-1])
s_saf_LN0= np.copy(wcerTS); s_saf_LN0[~((s_saf_ind==1)&(LN0ind==1))]=np.nan;s_saf_LN0=np.reshape(s_saf_LN0[6:-6],[12,numYrs-1])
s_saf_EN1= np.copy(wcerTS); s_saf_EN1[~((s_saf_ind==1)&(EN1ind==1))]=np.nan;s_saf_EN1=np.reshape(s_saf_EN1[6:-6],[12,numYrs-1])
s_saf_LN1= np.copy(wcerTS); s_saf_LN1[~((s_saf_ind==1)&(LN1ind==1))]=np.nan;s_saf_LN1=np.reshape(s_saf_LN1[6:-6],[12,numYrs-1])

s_saf_EN_1=s_saf_EN_1[~np.isnan(s_saf_EN_1)];s_saf_LN_1=s_saf_LN_1[~np.isnan(s_saf_LN_1)]
s_saf_EN0=s_saf_EN0[~np.isnan(s_saf_EN0)];s_saf_LN0=s_saf_LN0[~np.isnan(s_saf_LN0)]
s_saf_EN1=s_saf_EN1[~np.isnan(s_saf_EN1)];s_saf_LN1=s_saf_LN1[~np.isnan(s_saf_LN1)]
m_saf_EN_1=m_saf_EN_1[~np.isnan(m_saf_EN_1)];m_saf_LN_1=m_saf_LN_1[~np.isnan(m_saf_LN_1)]
m_saf_EN0=m_saf_EN0[~np.isnan(m_saf_EN0)];m_saf_LN0=m_saf_LN0[~np.isnan(m_saf_LN0)]
m_saf_EN1=m_saf_EN1[~np.isnan(m_saf_EN1)];m_saf_LN1=m_saf_LN1[~np.isnan(m_saf_LN1)]

s_EN_1 = np.append(np.append(s1_EN_1,s2_EN_1,0),s3_EN_1,0);s_EN_1=s_EN_1[~np.isnan(s_EN_1)]
s_LN_1 = np.append(np.append(s1_LN_1,s2_LN_1,0),s3_LN_1,0);s_LN_1=s_LN_1[~np.isnan(s_LN_1)]
s_EN0 = np.append(np.append(s1_EN0,s2_EN0,0),s3_EN0,0);s_EN0=s_EN0[~np.isnan(s_EN0)]
s_LN0 = np.append(np.append(s1_LN0,s2_LN0,0),s3_LN0,0);s_LN0=s_LN0[~np.isnan(s_LN0)]
s_EN1 = np.append(np.append(s1_EN1,s2_EN1,0),s3_EN1,0);s_EN1=s_EN1[~np.isnan(s_EN1)]
s_LN1 = np.append(np.append(s1_LN1,s2_LN1,0),s3_LN1,0);s_LN1=s_LN1[~np.isnan(s_LN1)]

m_EN_1 = np.append(np.append(m1_EN_1,m2_EN_1,0),m3_EN_1,0);m_EN_1=m_EN_1[~np.isnan(m_EN_1)]
m_LN_1 = np.append(np.append(m1_LN_1,m2_LN_1,0),m3_LN_1,0);m_LN_1=m_LN_1[~np.isnan(m_LN_1)]
m_EN0 = np.append(np.append(m1_EN0,m2_EN0,0),m3_EN0,0);m_EN0=m_EN0[~np.isnan(m_EN0)]
m_LN0 = np.append(np.append(m1_LN0,m2_LN0,0),m3_LN0,0);m_LN0=m_LN0[~np.isnan(m_LN0)]
m_EN1 = np.append(np.append(m1_EN1,m2_EN1,0),m3_EN1,0);m_EN1=m_EN1[~np.isnan(m_EN1)]
m_LN1 = np.append(np.append(m1_LN1,m2_LN1,0),m3_LN1,0);m_LN1=m_LN1[~np.isnan(m_LN1)]


# function for setting the colors of the box plots pairs
def setBoxColors(bp, boxNum, colorNam):
    setp(bp['boxes'][boxNum-1], color=colorNam)
    setp(bp['caps'][(boxNum*2)-2], color=colorNam)
    setp(bp['caps'][boxNum*2-1], color=colorNam)
    setp(bp['whiskers'][boxNum*2-2], color=colorNam)
    setp(bp['whiskers'][boxNum*2-1], color=colorNam)
    setp(bp['medians'][boxNum-1], color=colorNam)
#    setp(bp['fliers'][boxNum*2-2], color=colorNam)
#    setp(bp['fliers'][boxNum*2-1], color=colorNam)




fig = plt.figure(); ax1=plt.subplot(211);ax2=plt.subplot(212);
enBP = ax1.boxplot([s_EN_1,s_saf_EN_1, s_EN0,s_saf_EN0, s_EN1,s_saf_EN1],positions = [1,2,4,5,7,8],showmeans=True)
lnBP = ax2.boxplot([s_LN_1,s_saf_LN_1, s_LN0,s_saf_LN0, s_LN1,s_saf_LN1],positions = [1,2,4,5,7,8],showmeans=True)
setBoxColors(enBP, 1, 'darkblue');setBoxColors(enBP, 3, 'darkblue');setBoxColors(enBP, 5, 'darkblue')
setBoxColors(enBP, 2, 'red');setBoxColors(enBP, 4, 'red');setBoxColors(enBP, 6, 'red')
ax1.set_xticklabels(['EN -1','EN 0','EN 1']);ax2.set_xticklabels(['LN -1','LN 0','LN 1'])
ax1.set_xticks([1.5,4.5,7.5]);ax2.set_xticks([1.5,4.5,7.5])
setBoxColors(lnBP, 1, 'darkblue');setBoxColors(lnBP, 3, 'darkblue');setBoxColors(lnBP, 5, 'darkblue')
setBoxColors(lnBP, 2, 'red');setBoxColors(lnBP, 4, 'red');setBoxColors(lnBP, 6, 'red')
fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/soybean_SWC_'+seas+'boxPlot.png')
plt.close()


fig = plt.figure(); ax1=plt.subplot(211);ax2=plt.subplot(212);
enBP = ax1.boxplot([m_EN_1,m_saf_EN_1, m_EN0,m_saf_EN0, m_EN1,m_saf_EN1],positions = [1,2,4,5,7,8],showmeans=True)
lnBP = ax2.boxplot([m_LN_1,m_saf_LN_1, m_LN0,m_saf_LN0, m_LN1,m_saf_LN1],positions = [1,2,4,5,7,8],showmeans=True)
setBoxColors(enBP, 1, 'darkblue');setBoxColors(enBP, 3, 'darkblue');setBoxColors(enBP, 5, 'darkblue')
setBoxColors(enBP, 2, 'red');setBoxColors(enBP, 4, 'red');setBoxColors(enBP, 6, 'red')
ax1.set_xticklabels(['EN -1','EN 0','EN 1']);ax2.set_xticklabels(['LN -1','LN 0','LN 1'])
ax1.set_xticks([1.5,4.5,7.5]);ax2.set_xticks([1.5,4.5,7.5])
setBoxColors(lnBP, 1, 'darkblue');setBoxColors(lnBP, 3, 'darkblue');setBoxColors(lnBP, 5, 'darkblue')
setBoxColors(lnBP, 2, 'red');setBoxColors(lnBP, 4, 'red');setBoxColors(lnBP, 6, 'red')
fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/maize_SWC_'+seas+'boxPlot.png')
plt.close()