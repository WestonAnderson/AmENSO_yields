# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:24:43 2016

@author: westonanderson

Makes a call to the Soil water balance script, then plots the results as a 
spatial composite of EN and LN anomalies
"""
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from scipy.stats import ttest_ind
from scipy.stats import ttest_rel

years = [1950,2011] 
numYrs = years[1]-years[0]
ENyrs = np.array([1953, 1957, 1963, 1965, 1968, 1972, 1976, 1979, 1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009]) #EN
LNyrs = np.array([1954, 1964, 1970, 1973, 1983, 1988, 1995, 1998, 2007,2010]) #LN
seas = 'flwr'
plot = 'saf' #anom or absolute or saf or trd

ENind = np.repeat(np.array(ENyrs - years[0])*12,12) + np.tile(range(12),np.size(ENyrs))
ENyrInd = np.zeros([numYrs*12],dtype='int')
EN0ind =np.copy(ENyrInd); EN0ind[ENind] =1
EN_1ind =np.copy(ENyrInd); EN_1ind[ENind-12] =1
EN1ind =np.copy(ENyrInd); EN1ind[ENind+12] =1
EN2ind =np.copy(ENyrInd); EN2ind[ENind[ENind<708]+24] =1
#NOTE: the extra '+12' (+24) is because an EN calendar year ends in dec, but soybean/maize season starts Jan of next calendar year


LNind = np.repeat(np.array(LNyrs - years[0])*12,12) + np.tile(range(12),np.size(LNyrs))
LNyrInd = np.zeros([numYrs*12],dtype='int')
LN0ind =np.copy(LNyrInd); LN0ind[LNind] =1
LN_1ind =np.copy(LNyrInd); LN_1ind[LNind-12] =1
LN1ind =np.copy(LNyrInd); LN1ind[LNind[LNind<720]+12] =1
LN2ind =np.copy(LNyrInd); LN2ind[LNind[LNind<708]+24] =1

#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#
# Calculate the soil water balance
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/project/AmENSO_yields/soilWaterBalance.py')
#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#
lons,lats = np.meshgrid(lons,lats)
#Flowering seasons for each Kc curve by crop. Each season is two months long
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

#NOTE: the '+1' is because an EN calendar year ends in dec, 
#      but soybean/maize season starts jan or feb of next calendar year
s1_EN_1 = np.copy(w_s1);s1_EN_1[~((s1_ind==1)&(EN0ind==1)),...]=np.nan;s1_EN_1=np.reshape(s1_EN_1[6:-6,...],[12,numYrs-1,180,360])
s1_LN_1 = np.copy(w_s1);s1_LN_1[~((s1_ind==1)&(LN0ind==1)),...]=np.nan;s1_LN_1=np.reshape(s1_LN_1[6:-6,...],[12,numYrs-1,180,360])
s1_EN0 = np.copy(w_s1);s1_EN0[~((s1_ind==1)&(EN1ind==1)),...]=np.nan;s1_EN0=np.reshape(s1_EN0[6:-6,...],[12,numYrs-1,180,360])
s1_LN0 = np.copy(w_s1);s1_LN0[~((s1_ind==1)&(LN1ind==1)),...]=np.nan;s1_LN0=np.reshape(s1_LN0[6:-6,...],[12,numYrs-1,180,360])
s1_EN1 = np.copy(w_s1);s1_EN1[~((s1_ind==1)&(EN2ind==1)),...]=np.nan;s1_EN1=np.reshape(s1_EN1[6:-6,...],[12,numYrs-1,180,360])
s1_LN1 = np.copy(w_s1);s1_LN1[~((s1_ind==1)&(LN2ind==1)),...]=np.nan;s1_LN1=np.reshape(s1_LN1[6:-6,...],[12,numYrs-1,180,360])
s1_EN_1 = np.nanmean(s1_EN_1,0);s1_LN_1 = np.nanmean(s1_LN_1,0)
s1_EN0 = np.nanmean(s1_EN0,0);s1_LN0 = np.nanmean(s1_LN0,0)
s1_EN1 = np.nanmean(s1_EN1,0);s1_LN1 = np.nanmean(s1_LN1,0)

s2_EN_1 = np.copy(w_s2);s2_EN_1[~((s2_ind==1)&(EN0ind==1)),...]=np.nan;s2_EN_1=np.reshape(s2_EN_1[6:-6,...],[12,numYrs-1,180,360])
s2_LN_1 = np.copy(w_s2);s2_LN_1[~((s2_ind==1)&(LN0ind==1)),...]=np.nan;s2_LN_1=np.reshape(s2_LN_1[6:-6,...],[12,numYrs-1,180,360])
s2_EN0 = np.copy(w_s2);s2_EN0[~((s2_ind==1)&(EN1ind==1)),...]=np.nan;s2_EN0=np.reshape(s2_EN0[6:-6,...],[12,numYrs-1,180,360])
s2_LN0 = np.copy(w_s2);s2_LN0[~((s2_ind==1)&(LN1ind==1)),...]=np.nan;s2_LN0=np.reshape(s2_LN0[6:-6,...],[12,numYrs-1,180,360])
s2_EN1 = np.copy(w_s2);s2_EN1[~((s2_ind==1)&(EN2ind==1)),...]=np.nan;s2_EN1=np.reshape(s2_EN1[6:-6,...],[12,numYrs-1,180,360])
s2_LN1 = np.copy(w_s2);s2_LN1[~((s2_ind==1)&(LN2ind==1)),...]=np.nan;s2_LN1=np.reshape(s2_LN1[6:-6,...],[12,numYrs-1,180,360])
s2_EN_1 = np.nanmean(s2_EN_1,0);s2_LN_1 = np.nanmean(s2_LN_1,0)
s2_EN0 = np.nanmean(s2_EN0,0);s2_LN0 = np.nanmean(s2_LN0,0)
s2_EN1 = np.nanmean(s2_EN1,0);s2_LN1 = np.nanmean(s2_LN1,0)

s3_EN_1 = np.copy(w_s3);s3_EN_1[~((s3_ind==1)&(EN0ind==1)),...]=np.nan;s3_EN_1=np.reshape(s3_EN_1[6:-6,...],[12,numYrs-1,180,360])
s3_LN_1 = np.copy(w_s3);s3_LN_1[~((s3_ind==1)&(LN0ind==1)),...]=np.nan;s3_LN_1=np.reshape(s3_LN_1[6:-6,...],[12,numYrs-1,180,360])
s3_EN0 = np.copy(w_s3);s3_EN0[~((s3_ind==1)&(EN1ind==1)),...]=np.nan;s3_EN0=np.reshape(s3_EN0[6:-6,...],[12,numYrs-1,180,360])
s3_LN0 = np.copy(w_s3);s3_LN0[~((s3_ind==1)&(LN1ind==1)),...]=np.nan;s3_LN0=np.reshape(s3_LN0[6:-6,...],[12,numYrs-1,180,360])
s3_EN1 = np.copy(w_s3);s3_EN1[~((s3_ind==1)&(EN2ind==1)),...]=np.nan;s3_EN1=np.reshape(s3_EN1[6:-6,...],[12,numYrs-1,180,360])
s3_LN1 = np.copy(w_s3);s3_LN1[~((s3_ind==1)&(LN2ind==1)),...]=np.nan;s3_LN1=np.reshape(s3_LN1[6:-6,...],[12,numYrs-1,180,360])
s3_EN_1 = np.nanmean(s3_EN_1,0);s3_LN_1 = np.nanmean(s3_LN_1,0)
s3_EN0 = np.nanmean(s3_EN0,0);s3_LN0 = np.nanmean(s3_LN0,0)
s3_EN1 = np.nanmean(s3_EN1,0);s3_LN1 = np.nanmean(s3_LN1,0)

m1_EN_1 = np.copy(w_s1);m1_EN_1[~((m1_ind==1)&(EN0ind==1)),...]=np.nan;m1_EN_1=np.reshape(m1_EN_1[6:-6,...],[12,numYrs-1,180,360])
m1_LN_1 = np.copy(w_s1);m1_LN_1[~((m1_ind==1)&(LN0ind==1)),...]=np.nan;m1_LN_1=np.reshape(m1_LN_1[6:-6,...],[12,numYrs-1,180,360])
m1_EN0 = np.copy(w_s1);m1_EN0[~((m1_ind==1)&(EN1ind==1)),...]=np.nan;m1_EN0=np.reshape(m1_EN0[6:-6,...],[12,numYrs-1,180,360])
m1_LN0 = np.copy(w_s1);m1_LN0[~((m1_ind==1)&(LN1ind==1)),...]=np.nan;m1_LN0=np.reshape(m1_LN0[6:-6,...],[12,numYrs-1,180,360])
m1_EN1 = np.copy(w_s1);m1_EN1[~((m1_ind==1)&(EN2ind==1)),...]=np.nan;m1_EN1=np.reshape(m1_EN1[6:-6,...],[12,numYrs-1,180,360])
m1_LN1 = np.copy(w_s1);m1_LN1[~((m1_ind==1)&(LN2ind==1)),...]=np.nan;m1_LN1=np.reshape(m1_LN1[6:-6,...],[12,numYrs-1,180,360])
m1_EN_1 = np.nanmean(m1_EN_1,0);m1_LN_1 = np.nanmean(m1_LN_1,0)
m1_EN0 = np.nanmean(m1_EN0,0);m1_LN0 = np.nanmean(m1_LN0,0)
m1_EN1 = np.nanmean(m1_EN1,0);m1_LN1 = np.nanmean(m1_LN1,0)

m2_EN_1 = np.copy(w_s2);m2_EN_1[~((m2_ind==1)&(EN0ind==1)),...]=np.nan;m2_EN_1=np.reshape(m2_EN_1[6:-6,...],[12,numYrs-1,180,360])
m2_LN_1 = np.copy(w_s2);m2_LN_1[~((m2_ind==1)&(LN0ind==1)),...]=np.nan;m2_LN_1=np.reshape(m2_LN_1[6:-6,...],[12,numYrs-1,180,360])
m2_EN0 = np.copy(w_s2);m2_EN0[~((m2_ind==1)&(EN1ind==1)),...]=np.nan;m2_EN0=np.reshape(m2_EN0[6:-6,...],[12,numYrs-1,180,360])
m2_LN0 = np.copy(w_s2);m2_LN0[~((m2_ind==1)&(LN1ind==1)),...]=np.nan;m2_LN0=np.reshape(m2_LN0[6:-6,...],[12,numYrs-1,180,360])
m2_EN1 = np.copy(w_s2);m2_EN1[~((m2_ind==1)&(EN2ind==1)),...]=np.nan;m2_EN1=np.reshape(m2_EN1[6:-6,...],[12,numYrs-1,180,360])
m2_LN1 = np.copy(w_s2);m2_LN1[~((m2_ind==1)&(LN2ind==1)),...]=np.nan;m2_LN1=np.reshape(m2_LN1[6:-6,...],[12,numYrs-1,180,360])
m2_EN_1 = np.nanmean(m2_EN_1,0);m2_LN_1 = np.nanmean(m2_LN_1,0)
m2_EN0 = np.nanmean(m2_EN0,0);m2_LN0 = np.nanmean(m2_LN0,0)
m2_EN1 = np.nanmean(m2_EN1,0);m2_LN1 = np.nanmean(m2_LN1,0)

m3_EN_1 = np.copy(w_s3);m3_EN_1[~((m3_ind==1)&(EN0ind==1)),...]=np.nan;m3_EN_1=np.reshape(m3_EN_1[6:-6,...],[12,numYrs-1,180,360])
m3_LN_1 = np.copy(w_s3);m3_LN_1[~((m3_ind==1)&(LN0ind==1)),...]=np.nan;m3_LN_1=np.reshape(m3_LN_1[6:-6,...],[12,numYrs-1,180,360])
m3_EN0 = np.copy(w_s3);m3_EN0[~((m3_ind==1)&(EN1ind==1)),...]=np.nan;m3_EN0=np.reshape(m3_EN0[6:-6,...],[12,numYrs-1,180,360])
m3_LN0 = np.copy(w_s3);m3_LN0[~((m3_ind==1)&(LN1ind==1)),...]=np.nan;m3_LN0=np.reshape(m3_LN0[6:-6,...],[12,numYrs-1,180,360])
m3_EN1 = np.copy(w_s3);m3_EN1[~((m3_ind==1)&(EN2ind==1)),...]=np.nan;m3_EN1=np.reshape(m3_EN1[6:-6,...],[12,numYrs-1,180,360])
m3_LN1 = np.copy(w_s3);m3_LN1[~((m3_ind==1)&(LN2ind==1)),...]=np.nan;m3_LN1=np.reshape(m3_LN1[6:-6,...],[12,numYrs-1,180,360])
m3_EN_1 = np.nanmean(m3_EN_1,0);m3_LN_1 = np.nanmean(m3_LN_1,0)
m3_EN0 = np.nanmean(m3_EN0,0);m3_LN0 = np.nanmean(m3_LN0,0)
m3_EN1 = np.nanmean(m3_EN1,0);m3_LN1 = np.nanmean(m3_LN1,0)

m_saf_EN_1= np.copy(w_cer); m_saf_EN_1[~((m_saf_ind==1)&(EN0ind==1)),...]=np.nan;m_saf_EN_1=np.reshape(m_saf_EN_1[6:-6,...],[12,numYrs-1,180,360])
m_saf_LN_1= np.copy(w_cer); m_saf_LN_1[~((m_saf_ind==1)&(LN0ind==1)),...]=np.nan;m_saf_LN_1=np.reshape(m_saf_LN_1[6:-6,...],[12,numYrs-1,180,360])
m_saf_EN0= np.copy(w_cer); m_saf_EN0[~((m_saf_ind==1)&(EN1ind==1)),...]=np.nan;m_saf_EN0=np.reshape(m_saf_EN0[6:-6,...],[12,numYrs-1,180,360])
m_saf_LN0= np.copy(w_cer); m_saf_LN0[~((m_saf_ind==1)&(LN1ind==1)),...]=np.nan;m_saf_LN0=np.reshape(m_saf_LN0[6:-6,...],[12,numYrs-1,180,360])
m_saf_EN1= np.copy(w_cer); m_saf_EN1[~((m_saf_ind==1)&(EN2ind==1)),...]=np.nan;m_saf_EN1=np.reshape(m_saf_EN1[6:-6,...],[12,numYrs-1,180,360])
m_saf_LN1= np.copy(w_cer); m_saf_LN1[~((m_saf_ind==1)&(LN2ind==1)),...]=np.nan;m_saf_LN1=np.reshape(m_saf_LN1[6:-6,...],[12,numYrs-1,180,360])
m_saf_EN_1 = np.nanmean(m_saf_EN_1,0);m_saf_LN_1 = np.nanmean(m_saf_LN_1,0)
m_saf_EN0 = np.nanmean(m_saf_EN0,0);m_saf_LN0 = np.nanmean(m_saf_LN0,0)
m_saf_EN1 = np.nanmean(m_saf_EN1,0);m_saf_LN1 = np.nanmean(m_saf_LN1,0)

#Safrinha soy is the only anomaly that happens in the same 'calendar year' as the enso event
s_saf_EN_1= np.copy(w_cer); s_saf_EN_1[~((s_saf_ind==1)&(EN_1ind==1)),...]=np.nan;s_saf_EN_1=np.reshape(s_saf_EN_1[6:-6,...],[12,numYrs-1,180,360])
s_saf_LN_1= np.copy(w_cer); s_saf_LN_1[~((s_saf_ind==1)&(LN_1ind==1)),...]=np.nan;s_saf_LN_1=np.reshape(s_saf_LN_1[6:-6,...],[12,numYrs-1,180,360])
s_saf_EN0= np.copy(w_cer); s_saf_EN0[~((s_saf_ind==1)&(EN0ind==1)),...]=np.nan;s_saf_EN0=np.reshape(s_saf_EN0[6:-6,...],[12,numYrs-1,180,360])
s_saf_LN0= np.copy(w_cer); s_saf_LN0[~((s_saf_ind==1)&(LN0ind==1)),...]=np.nan;s_saf_LN0=np.reshape(s_saf_LN0[6:-6,...],[12,numYrs-1,180,360])
s_saf_EN1= np.copy(w_cer); s_saf_EN1[~((s_saf_ind==1)&(EN1ind==1)),...]=np.nan;s_saf_EN1=np.reshape(s_saf_EN1[6:-6,...],[12,numYrs-1,180,360])
s_saf_LN1= np.copy(w_cer); s_saf_LN1[~((s_saf_ind==1)&(LN1ind==1)),...]=np.nan;s_saf_LN1=np.reshape(s_saf_LN1[6:-6,...],[12,numYrs-1,180,360])
s_saf_EN_1 = np.nanmean(s_saf_EN_1,0);s_saf_LN_1 = np.nanmean(s_saf_LN_1,0)
s_saf_EN0 = np.nanmean(s_saf_EN0,0);s_saf_LN0 = np.nanmean(s_saf_LN0,0)
s_saf_EN1 = np.nanmean(s_saf_EN1,0);s_saf_LN1 = np.nanmean(s_saf_LN1,0)

s_EN_1 = np.append(np.append(s1_EN_1,s2_EN_1,0),s3_EN_1,0);
s_LN_1 = np.append(np.append(s1_LN_1,s2_LN_1,0),s3_LN_1,0);
s_EN0 = np.append(np.append(s1_EN0,s2_EN0,0),s3_EN0,0);
s_LN0 = np.append(np.append(s1_LN0,s2_LN0,0),s3_LN0,0);
s_EN1 = np.append(np.append(s1_EN1,s2_EN1,0),s3_EN1,0)
s_LN1 = np.append(np.append(s1_LN1,s2_LN1,0),s3_LN1,0)

m_EN_1 = np.append(np.append(m1_EN_1,m2_EN_1,0),m3_EN_1,0);
m_LN_1 = np.append(np.append(m1_LN_1,m2_LN_1,0),m3_LN_1,0);
m_EN0 = np.append(np.append(m1_EN0,m2_EN0,0),m3_EN0,0);
m_LN0 = np.append(np.append(m1_LN0,m2_LN0,0),m3_LN0,0);
m_EN1 = np.append(np.append(m1_EN1,m2_EN1,0),m3_EN1,0);
m_LN1 = np.append(np.append(m1_LN1,m2_LN1,0),m3_LN1,0);

#seasonal mean values
s1_tot = w_s1[(s1_ind==1),...];s2_tot = w_s2[(s2_ind==1),...];s3_tot = w_s3[(s3_ind==1),...];
m1_tot = w_s1[(m1_ind==1),...];m2_tot = w_s2[(m2_ind==1),...];m3_tot = w_s3[(m3_ind==1),...];
tot_s =  np.append(np.append(s1_tot,s2_tot,0),s3_tot,0)
tot_m =  np.append(np.append(m1_tot,m2_tot,0),m3_tot,0)

sSigEN_1 = np.zeros([180,360]);sSigLN_1 = np.zeros([180,360])
sSigEN0 = np.zeros([180,360]);sSigLN0 = np.zeros([180,360])
sSigEN1 = np.zeros([180,360]);sSigLN1 = np.zeros([180,360])
mSigEN_1 = np.zeros([180,360]);mSigLN_1 = np.zeros([180,360])
mSigEN0 = np.zeros([180,360]);mSigLN0 = np.zeros([180,360])
mSigEN1 = np.zeros([180,360]);mSigLN1 = np.zeros([180,360])
sSigTot = np.zeros([180,360]);mSigTot = np.zeros([180,360])

tot_s = np.nanmean(tot_s,0);tot_m = np.nanmean(tot_m,0)
saf_tot_s = np.nanmean(w_cer[(s_saf_ind==1),...],0)
saf_tot_m = np.nanmean(w_cer[(m_saf_ind==1),...],0)

plt.ioff()
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#define a helper function to pull up to date boundaries for Brazil
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_boundaries.py')
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~
######## Make maps of the SWC anomalies ###############

ur_lat2 = -3;ll_lat2 = -40;ll_lon2 = -70.;ur_lon2 = -35.
ur_lat = -3.;ll_lat = -40;ll_lon = -70.;ur_lon = -35.




if plot is 'saf':
    #test for significance in south america
    sSigEN_1 = ttest_ind(s_EN_1-tot_s,s_saf_EN_1-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN_1 = ttest_ind(s_LN_1-tot_s,s_saf_LN_1-saf_tot_s,nan_policy='omit').pvalue
    sSigEN0 = ttest_ind(s_EN0-tot_s,s_saf_EN0-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN0 = ttest_ind(s_LN0-tot_s,s_saf_LN0-saf_tot_s,nan_policy='omit').pvalue
    sSigEN1 = ttest_ind(s_EN1-tot_s,s_saf_EN1-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN1 = ttest_ind(s_LN1-tot_s,s_saf_LN1-saf_tot_s,nan_policy='omit').pvalue
    mSigEN_1 = ttest_ind(m_EN_1-tot_m,m_saf_EN_1-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN_1 = ttest_ind(m_LN_1-tot_m,m_saf_LN_1-saf_tot_m,nan_policy='omit').pvalue
    mSigEN0 = ttest_ind(m_EN0-tot_m,m_saf_EN0-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN0 = ttest_ind(m_LN0-tot_m,m_saf_LN0-saf_tot_m,nan_policy='omit').pvalue
    mSigEN1 = ttest_ind(m_EN1-tot_m,m_saf_EN1-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN1 = ttest_ind(m_LN1-tot_m,m_saf_LN1-saf_tot_m,nan_policy='omit').pvalue
    sSigTot = ttest_ind(tot_s,w_cer[(s_saf_ind==1),...],nan_policy='omit').pvalue
    mSigTot = ttest_ind(tot_m,w_cer[(m_saf_ind==1),...],nan_policy='omit').pvalue
    
    
    s_EN_1 = np.nanmean(s_EN_1,0);s_LN_1 = np.nanmean(s_LN_1,0)
    m_EN_1 = np.nanmean(m_EN_1,0);m_LN_1 = np.nanmean(m_LN_1,0)
    s_saf_EN_1 = np.nanmean(s_saf_EN_1,0);s_saf_LN_1 = np.nanmean(s_saf_LN_1,0)
    m_saf_EN_1 = np.nanmean(m_saf_EN_1,0);m_saf_LN_1 = np.nanmean(m_saf_LN_1,0)
    s_EN0 = np.nanmean(s_EN0,0);s_LN0 = np.nanmean(s_LN0,0)
    m_EN0 = np.nanmean(m_EN0,0);m_LN0 = np.nanmean(m_LN0,0)
    s_saf_EN0 = np.nanmean(s_saf_EN0,0);s_saf_LN0 = np.nanmean(s_saf_LN0,0)
    m_saf_EN0 = np.nanmean(m_saf_EN0,0);m_saf_LN0 = np.nanmean(m_saf_LN0,0)
    s_EN1 = np.nanmean(s_EN1,0);s_LN1 = np.nanmean(s_LN1,0)
    m_EN1 = np.nanmean(m_EN1,0);m_LN1 = np.nanmean(m_LN1,0)
    s_saf_EN1 = np.nanmean(s_saf_EN1,0);s_saf_LN1 = np.nanmean(s_saf_LN1,0)
    m_saf_EN1 = np.nanmean(m_saf_EN1,0);m_saf_LN1 = np.nanmean(m_saf_LN1,0)

    totDiff_m = saf_tot_m-tot_m; totDiff_s = saf_tot_s-tot_s; 
    totDiff_m[totDiff_m>=.20]=.20;totDiff_m[totDiff_m<=-.20]=-.20
    totDiff_s[totDiff_s>=.20]=.20;totDiff_s[totDiff_s<=-.20]=-.20
    clevs = np.arange(-.20,.21,.01)

    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC'_{saf}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC'_{saf}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,totDiff_m,clevs,shading='flat',cmap='PuOr',latlon=True)
    m1.drawcoastlines();m1.drawcountries()
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(m_saf_EN_1-saf_tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(m_saf_EN0-saf_tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(m_saf_EN1-saf_tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(m_saf_LN_1-saf_tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(m_saf_LN0-saf_tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(m_saf_LN1-saf_tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m13.drawcoastlines();m13.drawcountries();
    
    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/maize_SWC_'+seas+'Saf.png')
    
    
    
    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC'_{saf}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC'_{saf}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,totDiff_s,clevs,shading='flat',cmap='PuOr',latlon=True)
    m1.drawcoastlines();m1.drawcountries()
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(s_saf_EN_1-saf_tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(s_saf_EN0-saf_tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(s_saf_EN1-saf_tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(s_saf_LN_1-saf_tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(s_saf_LN0-saf_tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(s_saf_LN1-saf_tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m13.drawcoastlines();m13.drawcountries();

    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/soybean_SWC_'+seas+'Saf.png')


if plot is 'trd':
    #test for significance in south america
    sSigEN_1 = ttest_ind(s_EN_1-tot_s,s_saf_EN_1-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN_1 = ttest_ind(s_LN_1-tot_s,s_saf_LN_1-saf_tot_s,nan_policy='omit').pvalue
    sSigEN0 = ttest_ind(s_EN0-tot_s,s_saf_EN0-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN0 = ttest_ind(s_LN0-tot_s,s_saf_LN0-saf_tot_s,nan_policy='omit').pvalue
    sSigEN1 = ttest_ind(s_EN1-tot_s,s_saf_EN1-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN1 = ttest_ind(s_LN1-tot_s,s_saf_LN1-saf_tot_s,nan_policy='omit').pvalue
    mSigEN_1 = ttest_ind(m_EN_1-tot_m,m_saf_EN_1-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN_1 = ttest_ind(m_LN_1-tot_m,m_saf_LN_1-saf_tot_m,nan_policy='omit').pvalue
    mSigEN0 = ttest_ind(m_EN0-tot_m,m_saf_EN0-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN0 = ttest_ind(m_LN0-tot_m,m_saf_LN0-saf_tot_m,nan_policy='omit').pvalue
    mSigEN1 = ttest_ind(m_EN1-tot_m,m_saf_EN1-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN1 = ttest_ind(m_LN1-tot_m,m_saf_LN1-saf_tot_m,nan_policy='omit').pvalue
    sSigTot = ttest_ind(tot_s,w_cer[(s_saf_ind==1),...],nan_policy='omit').pvalue
    mSigTot = ttest_ind(tot_m,w_cer[(m_saf_ind==1),...],nan_policy='omit').pvalue
    
    
    s_EN_1 = np.nanmean(s_EN_1,0);s_LN_1 = np.nanmean(s_LN_1,0)
    m_EN_1 = np.nanmean(m_EN_1,0);m_LN_1 = np.nanmean(m_LN_1,0)
    s_saf_EN_1 = np.nanmean(s_saf_EN_1,0);s_saf_LN_1 = np.nanmean(s_saf_LN_1,0)
    m_saf_EN_1 = np.nanmean(m_saf_EN_1,0);m_saf_LN_1 = np.nanmean(m_saf_LN_1,0)
    s_EN0 = np.nanmean(s_EN0,0);s_LN0 = np.nanmean(s_LN0,0)
    m_EN0 = np.nanmean(m_EN0,0);m_LN0 = np.nanmean(m_LN0,0)
    s_saf_EN0 = np.nanmean(s_saf_EN0,0);s_saf_LN0 = np.nanmean(s_saf_LN0,0)
    m_saf_EN0 = np.nanmean(m_saf_EN0,0);m_saf_LN0 = np.nanmean(m_saf_LN0,0)
    s_EN1 = np.nanmean(s_EN1,0);s_LN1 = np.nanmean(s_LN1,0)
    m_EN1 = np.nanmean(m_EN1,0);m_LN1 = np.nanmean(m_LN1,0)
    s_saf_EN1 = np.nanmean(s_saf_EN1,0);s_saf_LN1 = np.nanmean(s_saf_LN1,0)
    m_saf_EN1 = np.nanmean(m_saf_EN1,0);m_saf_LN1 = np.nanmean(m_saf_LN1,0)


    clevs = np.arange(-.20,.21,.01)

    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC'_{trd}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC'_{trd}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,saf_tot_m-tot_m,clevs,shading='flat',cmap='PuOr',latlon=True)
    m1.drawcoastlines();m1.drawcountries()
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(m_EN_1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(m_EN0-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(m_EN1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(m_LN_1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(m_LN0-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(m_LN1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    m13.drawcoastlines();m13.drawcountries();
    
    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/maize_SWC_'+seas+'Trd.png')
    
    
    
    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC'_{trd}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC'_{trd}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,saf_tot_s-tot_s,clevs,shading='flat',cmap='PuOr',latlon=True)
    m1.drawcoastlines();m1.drawcountries()
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(s_EN_1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(s_EN0-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(s_EN1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(s_LN_1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(s_LN0-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(s_LN1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    m13.drawcoastlines();m13.drawcountries();

    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/soybean_SWC_'+seas+'Trd.png')



if plot is 'anom':
    #test for significance in south america
    sSigEN_1 = ttest_ind(s_EN_1-tot_s,s_saf_EN_1-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN_1 = ttest_ind(s_LN_1-tot_s,s_saf_LN_1-saf_tot_s,nan_policy='omit').pvalue
    sSigEN0 = ttest_ind(s_EN0-tot_s,s_saf_EN0-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN0 = ttest_ind(s_LN0-tot_s,s_saf_LN0-saf_tot_s,nan_policy='omit').pvalue
    sSigEN1 = ttest_ind(s_EN1-tot_s,s_saf_EN1-saf_tot_s,nan_policy='omit').pvalue ;
    sSigLN1 = ttest_ind(s_LN1-tot_s,s_saf_LN1-saf_tot_s,nan_policy='omit').pvalue
    mSigEN_1 = ttest_ind(m_EN_1-tot_m,m_saf_EN_1-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN_1 = ttest_ind(m_LN_1-tot_m,m_saf_LN_1-saf_tot_m,nan_policy='omit').pvalue
    mSigEN0 = ttest_ind(m_EN0-tot_m,m_saf_EN0-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN0 = ttest_ind(m_LN0-tot_m,m_saf_LN0-saf_tot_m,nan_policy='omit').pvalue
    mSigEN1 = ttest_ind(m_EN1-tot_m,m_saf_EN1-saf_tot_m,nan_policy='omit').pvalue ;
    mSigLN1 = ttest_ind(m_LN1-tot_m,m_saf_LN1-saf_tot_m,nan_policy='omit').pvalue
    sSigTot = ttest_ind(tot_s,w_cer[(s_saf_ind==1),...],nan_policy='omit').pvalue
    mSigTot = ttest_ind(tot_m,w_cer[(m_saf_ind==1),...],nan_policy='omit').pvalue
    
    
    s_EN_1 = np.nanmean(s_EN_1,0);s_LN_1 = np.nanmean(s_LN_1,0)
    m_EN_1 = np.nanmean(m_EN_1,0);m_LN_1 = np.nanmean(m_LN_1,0)
    s_saf_EN_1 = np.nanmean(s_saf_EN_1,0);s_saf_LN_1 = np.nanmean(s_saf_LN_1,0)
    m_saf_EN_1 = np.nanmean(m_saf_EN_1,0);m_saf_LN_1 = np.nanmean(m_saf_LN_1,0)
    s_EN0 = np.nanmean(s_EN0,0);s_LN0 = np.nanmean(s_LN0,0)
    m_EN0 = np.nanmean(m_EN0,0);m_LN0 = np.nanmean(m_LN0,0)
    s_saf_EN0 = np.nanmean(s_saf_EN0,0);s_saf_LN0 = np.nanmean(s_saf_LN0,0)
    m_saf_EN0 = np.nanmean(m_saf_EN0,0);m_saf_LN0 = np.nanmean(m_saf_LN0,0)
    s_EN1 = np.nanmean(s_EN1,0);s_LN1 = np.nanmean(s_LN1,0)
    m_EN1 = np.nanmean(m_EN1,0);m_LN1 = np.nanmean(m_LN1,0)
    s_saf_EN1 = np.nanmean(s_saf_EN1,0);s_saf_LN1 = np.nanmean(s_saf_LN1,0)
    m_saf_EN1 = np.nanmean(m_saf_EN1,0);m_saf_LN1 = np.nanmean(m_saf_LN1,0)


    clevs = np.arange(-.20,.21,.01)

    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC'_{saf}-SWC'_{trd}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC'_{saf}-SWC'_{trd}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,saf_tot_m-tot_m,clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m1(lons[mSigTot<=0.1],lats[mSigTot<=0.1]);m1.plot(sigLons,sigLats,'ko',markersize=0.5)
    m1.drawcoastlines();m1.drawcountries()
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(m_saf_EN_1-saf_tot_m)-(m_EN_1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m21(lons[mSigEN_1<=0.1],lats[mSigEN_1<=0.1]);m21.plot(sigLons,sigLats,'ko',markersize=0.5)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(m_saf_EN0-saf_tot_m)-(m_EN0-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m22(lons[mSigEN0<=0.1],lats[mSigEN0<=0.1]);m22.plot(sigLons,sigLats,'ko',markersize=0.5)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(m_saf_EN1-saf_tot_m)-(m_EN1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m23(lons[mSigEN1<=0.1],lats[mSigEN1<=0.1]);m23.plot(sigLons,sigLats,'ko',markersize=0.5)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(m_saf_LN_1-saf_tot_m)-(m_LN_1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m11(lons[mSigLN_1<=0.1],lats[mSigLN_1<=0.1]);m11.plot(sigLons,sigLats,'ko',markersize=0.5)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(m_saf_LN0-saf_tot_m)-(m_LN0-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m12(lons[mSigLN0<=0.1],lats[mSigLN0<=0.1]);m12.plot(sigLons,sigLats,'ko',markersize=0.5)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(m_saf_LN1-saf_tot_m)-(m_LN1-tot_m),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m13(lons[mSigLN1<=0.1],lats[mSigLN1<=0.1]);m13.plot(sigLons,sigLats,'ko',markersize=0.5)
    m13.drawcoastlines();m13.drawcountries();
    
    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/maize_SWC_'+seas+'Anom.png')
    
    
    
    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC'_{saf}-SWC'_{trd}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC'_{saf}-SWC'_{trd}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,saf_tot_s-tot_s,clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m1(lons[sSigTot<=0.1],lats[sSigTot<=0.1]);m1.plot(sigLons,sigLats,'ko',markersize=0.5)
    m1.drawcoastlines();m1.drawcountries()
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(s_saf_EN_1-saf_tot_s)-(s_EN_1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m21(lons[sSigEN_1<=0.1],lats[sSigEN_1<=0.1]);m21.plot(sigLons,sigLats,'ko',markersize=0.5)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(s_saf_EN0-saf_tot_s)-(s_EN0-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m22(lons[sSigEN0<=0.1],lats[sSigEN0<=0.1]);m22.plot(sigLons,sigLats,'ko',markersize=0.5)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(s_saf_EN1-saf_tot_s)-(s_EN1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m23(lons[sSigEN1<=0.1],lats[sSigEN1<=0.1]);m23.plot(sigLons,sigLats,'ko',markersize=0.5)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(s_saf_LN_1-saf_tot_s)-(s_LN_1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m11(lons[sSigLN_1<=0.1],lats[sSigLN_1<=0.1]);m11.plot(sigLons,sigLats,'ko',markersize=0.5)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(s_saf_LN0-saf_tot_s)-(s_LN0-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m12(lons[sSigLN0<=0.1],lats[sSigLN0<=0.1]);m12.plot(sigLons,sigLats,'ko',markersize=0.5)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(s_saf_LN1-saf_tot_s)-(s_LN1-tot_s),clevs,shading='flat',cmap='PuOr',latlon=True)
    sigLons, sigLats = m13(lons[sSigLN1<=0.1],lats[sSigLN1<=0.1]);m13.plot(sigLons,sigLats,'ko',markersize=0.5)
    m13.drawcoastlines();m13.drawcountries();

    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/soybean_SWC_'+seas+'Anom.png')

elif plot is 'absolute':
    #test for significance in south america
    sSigEN_1 = ttest_ind(s_EN_1,s_saf_EN_1,nan_policy='omit').pvalue ;
    sSigLN_1 = ttest_ind(s_LN_1,s_saf_LN_1,nan_policy='omit').pvalue
    sSigEN0 = ttest_ind(s_EN0,s_saf_EN0,nan_policy='omit').pvalue ;
    sSigLN0 = ttest_ind(s_LN0,s_saf_LN0,nan_policy='omit').pvalue
    sSigEN1 = ttest_ind(s_EN1,s_saf_EN1,nan_policy='omit').pvalue ;
    sSigLN1 = ttest_ind(s_LN1,s_saf_LN1,nan_policy='omit').pvalue
    mSigEN_1 = ttest_ind(m_EN_1,m_saf_EN_1,nan_policy='omit').pvalue ;
    mSigLN_1 = ttest_ind(m_LN_1,m_saf_LN_1,nan_policy='omit').pvalue
    mSigEN0 = ttest_ind(m_EN0,m_saf_EN0,nan_policy='omit').pvalue ;
    mSigLN0 = ttest_ind(m_LN0,m_saf_LN0,nan_policy='omit').pvalue
    mSigEN1 = ttest_ind(m_EN1,m_saf_EN1,nan_policy='omit').pvalue ;
    mSigLN1 = ttest_ind(m_LN1,m_saf_LN1,nan_policy='omit').pvalue
    sSigTot = ttest_ind(tot_s,w_cer[(s_saf_ind==1),...],nan_policy='omit').pvalue
    mSigTot = ttest_ind(tot_m,w_cer[(m_saf_ind==1),...],nan_policy='omit').pvalue
    
    
    s_EN_1 = np.nanmean(s_EN_1,0);s_LN_1 = np.nanmean(s_LN_1,0)
    m_EN_1 = np.nanmean(m_EN_1,0);m_LN_1 = np.nanmean(m_LN_1,0)
    s_saf_EN_1 = np.nanmean(s_saf_EN_1,0);s_saf_LN_1 = np.nanmean(s_saf_LN_1,0)
    m_saf_EN_1 = np.nanmean(m_saf_EN_1,0);m_saf_LN_1 = np.nanmean(m_saf_LN_1,0)
    s_EN0 = np.nanmean(s_EN0,0);s_LN0 = np.nanmean(s_LN0,0)
    m_EN0 = np.nanmean(m_EN0,0);m_LN0 = np.nanmean(m_LN0,0)
    s_saf_EN0 = np.nanmean(s_saf_EN0,0);s_saf_LN0 = np.nanmean(s_saf_LN0,0)
    m_saf_EN0 = np.nanmean(m_saf_EN0,0);m_saf_LN0 = np.nanmean(m_saf_LN0,0)
    s_EN1 = np.nanmean(s_EN1,0);s_LN1 = np.nanmean(s_LN1,0)
    m_EN1 = np.nanmean(m_EN1,0);m_LN1 = np.nanmean(m_LN1,0)
    s_saf_EN1 = np.nanmean(s_saf_EN1,0);s_saf_LN1 = np.nanmean(s_saf_LN1,0)
    m_saf_EN1 = np.nanmean(m_saf_EN1,0);m_saf_LN1 = np.nanmean(m_saf_LN1,0)
    
        
    clevs = np.arange(-.40,.41,.01)
    
    
    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC_{saf,EN}-SWC_{trd,EN}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC_{saf,LN}-SWC_{trd,LN}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,saf_tot_m-tot_m,clevs,shading='flat',cmap='PuOr',latlon=True)
    m1.drawcoastlines();m1.drawcountries()
    
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(m_saf_EN_1)-(m_EN_1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(m_saf_EN0)-(m_EN0),clevs,shading='flat',cmap='PuOr',latlon=True)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(m_saf_EN1)-(m_EN1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(m_saf_LN_1)-(m_LN_1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(m_saf_LN0)-(m_LN0),clevs,shading='flat',cmap='PuOr',latlon=True)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(m_saf_LN1)-(m_LN1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m13.drawcoastlines();m13.drawcountries();

    sigLons, sigLats = m1(lons[mSigTot<=0.1],lats[mSigTot<=0.1]);m1.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m21(lons[mSigEN_1<=0.1],lats[mSigEN_1<=0.1]);m21.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m22(lons[mSigEN0<=0.1],lats[mSigEN0<=0.1]);m22.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m23(lons[mSigEN1<=0.1],lats[mSigEN1<=0.1]);m23.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m11(lons[mSigLN_1<=0.1],lats[mSigLN_1<=0.1]);m11.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m12(lons[mSigLN0<=0.1],lats[mSigLN0<=0.1]);m12.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m13(lons[mSigLN1<=0.1],lats[mSigLN1<=0.1]);m13.plot(sigLons,sigLats,'ko',markersize=0.5)
    
    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/maize_SWC_'+seas+'.png')
    
    
    fig = plt.figure();
    #fig.suptitle(crp,fontsize=20)
    ax1=plt.subplot(221);
    ax21=plt.subplot(343);ax22=plt.subplot(347);ax23=plt.subplot(3,4,11);
    ax11=plt.subplot(344);ax12=plt.subplot(348);ax13=plt.subplot(3,4,12);
    ax1.set_title(u"\n $\overline{SWC}_{saf}-\overline{SWC}_{trd}$")
    ax11.set_title(u"La Niña life-cycle\n $SWC_{saf,LN}-SWC_{trd,LN}$", fontsize = 10)
    ax21.set_title(u"El Niño life-cycle\n $SWC_{saf,EN}-SWC_{trd,EN}$", fontsize = 10)
    ax11.yaxis.labelpad = 15;ax12.yaxis.labelpad = 15;ax13.yaxis.labelpad = 15
    y1 = ax11.set_ylabel('Year -1');y1.set_rotation(-90);ax11.yaxis.set_label_position('right')
    y2 = ax12.set_ylabel('Year 0');y2.set_rotation(-90);ax12.yaxis.set_label_position('right')
    y3 = ax13.set_ylabel('Year +1');y3.set_rotation(-90);ax13.yaxis.set_label_position('right')
    
    m1 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1);
    im1=m1.contourf(lons,lats,saf_tot_s-tot_s,clevs,shading='flat',cmap='PuOr',latlon=True)
    m1.drawcoastlines();m1.drawcountries()
    
    m21 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax21);
    m21.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im21=m21.contourf(lons,lats,(s_saf_EN_1)-(s_EN_1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax22);
    m22.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im22=m22.contourf(lons,lats,(s_saf_EN0)-(s_EN0),clevs,shading='flat',cmap='PuOr',latlon=True)
    m22.drawcoastlines();m22.drawcountries()
    m23 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax23);
    m23.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im23=m23.contourf(lons,lats,(s_saf_EN1)-(s_EN1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m23.drawcoastlines();m23.drawcountries()
    
    m11 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax11);
    m11.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im11=m11.contourf(lons,lats,(s_saf_LN_1)-(s_LN_1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m11.drawcoastlines();m11.drawcountries();
    m12 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax12);
    m12.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im12=m12.contourf(lons,lats,(s_saf_LN0)-(s_LN0),clevs,shading='flat',cmap='PuOr',latlon=True)
    m12.drawcoastlines();m12.drawcountries();
    m13 = Basemap(projection='mill',lon_0=lon0,lat_0=0,llcrnrlat=ll_lat,urcrnrlat=ur_lat,llcrnrlon=ll_lon,urcrnrlon=ur_lon,ax=ax13);
    m13.drawlsmask(resolution='c',land_color='0.8', ocean_color='w')
    im13=m13.contourf(lons,lats,(s_saf_LN1)-(s_LN1),clevs,shading='flat',cmap='PuOr',latlon=True)
    m13.drawcoastlines();m13.drawcountries();

    sigLons, sigLats = m1(lons[sSigTot<=0.1],lats[sSigTot<=0.1]);m1.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m21(lons[sSigEN_1<=0.1],lats[sSigEN_1<=0.1]);m21.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m22(lons[sSigEN0<=0.1],lats[sSigEN0<=0.1]);m22.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m23(lons[sSigEN1<=0.1],lats[sSigEN1<=0.1]);m23.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m11(lons[sSigLN_1<=0.1],lats[sSigLN_1<=0.1]);m11.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m12(lons[sSigLN0<=0.1],lats[sSigLN0<=0.1]);m12.plot(sigLons,sigLats,'ko',markersize=0.5)
    sigLons, sigLats = m13(lons[sSigLN1<=0.1],lats[sSigLN1<=0.1]);m13.plot(sigLons,sigLats,'ko',markersize=0.5)
    
    drawstates(m1,ax1,'BRA',1);
    drawstates(m21,ax21,'BRA',1);drawstates(m22,ax22,'BRA',1);drawstates(m23,ax23,'BRA',1)
    drawstates(m11,ax11,'BRA',1);drawstates(m12,ax12,'BRA',1);drawstates(m13,ax13,'BRA',1)
    
    cb1 = m1.colorbar(im21,"left", size="5%", pad="2%")
    cb1.ax.yaxis.set_ticks_position('left')
    ax1.legend(bbox_to_anchor=(.88, 1), loc=2, borderaxespad=2.)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/ENSOlifecycle/soybean_SWC_'+seas+'.png')
    

 