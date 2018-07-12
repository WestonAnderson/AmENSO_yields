# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:31:56 2016

@author: westonanderson

Makes a call to the Soil water balance script, then plots the results as a 
time series of EN and LN composites
"""
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon

plot = 'P_tot'
reg = 'SulCerr' #Sul, Cerr, GO, MS, RS or MT
crp = 'soybean'
ENyrs = np.array([1953, 1957, 1963, 1965, 1968, 1972, 1976, 1979, 1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009]) #EN
LNyrs = np.array([1954, 1964, 1970, 1973, 1983, 1988, 1995, 1998, 2007, 2010]) #LN


#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#
# Calculate the soil water balance
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/project/AmENSO_yields/soilWaterBalance.py')
#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#


#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~
######## Set the region of analysis ###############
if reg is 'SE_SA':
    #south region of Brazil
    lonMin = -60; lonMax = -52
    latMin = -21; latMax = -28
elif reg is 'Sul':
    #south region of Brazil
    lonMin = -57; lonMax = -47
    latMin = -22; latMax = -28
elif reg is 'Cerr':
    #West/Central Region
    lonMin = -55; lonMax = -45
    latMin =-15 ; latMax = -22
elif reg is 'RS':
    #State RS (early maize)
    lonMin = -55; lonMax = -50
    latMin = -27; latMax = -33
elif reg is 'MT':
    # state MT (safrinha)
    lonMin = -57; lonMax = -55
    latMin = -11; latMax = -14
elif reg is 'GO':
    # state GO (safrinha)
    lonMin = -53; lonMax = -49
    latMin = -16; latMax = -18
elif reg is 'MS':
    # state MS (safrinha)
    lonMin = -56.5; lonMax = -54.5
    latMin = -20; latMax = -23
elif reg is 'PR':
    lonMin = -54; lonMax = -48
    latMin = -22; latMax = -25
elif reg is 'MS_PR':
    lonMin = -57; lonMax = -47
    latMin = -20; latMax = -25
elif reg is 'SulCerr':
    lonMin = -57; lonMax = -47
    latMin = -15; latMax = -25    
    
latMax = np.min(np.where(lats<latMax))
latMin = np.max(np.where(lats>latMin))
lonMin = np.max(np.where(lons<lonMin))
lonMax = np.min(np.where(lons>lonMax))



#define ENSO years
EN_index = (ENyrs-years[0]-1)*12; LN_index = (LNyrs-years[0]-1)*12 # -1 b/c need to start the year before the ENSO year
EN_index = np.repeat(EN_index,42)+np.tile(range(42),EN_index.shape[0])
LN_index = np.repeat(LN_index,42)+np.tile(range(42),LN_index.shape[0])
#data doesn't go through 2011, so cut it off and fill with a blank year
LN_index = LN_index[LN_index<732]
LN_blank = np.zeros([18,180,360])*np.nan
EN_index = EN_index[EN_index<732]
EN_blank = np.zeros([6,180,360])*np.nan

#create composite objects
full_s1 = np.reshape(w_s1,[w_s1.shape[0]/12,12,w_s1.shape[1],w_s1.shape[2]]);
full_s2 = np.reshape(w_s2,[w_s2.shape[0]/12,12,w_s2.shape[1],w_s2.shape[2]]);
full_s3 = np.reshape(w_s3,[w_s3.shape[0]/12,12,w_s3.shape[1],w_s3.shape[2]]);
full_cer = np.reshape(w_cer,[w_cer.shape[0]/12,12,w_cer.shape[1],w_cer.shape[2]]);

full_P = np.reshape(p,[p.shape[0]/12,12,p.shape[1],p.shape[2]]); 
full_ETcer = np.reshape(Ecer,[Ecer.shape[0]/12,12,Ecer.shape[1],Ecer.shape[2]]); #full_ETs = np.zeros([12,180,360]); 
full_ET0 = np.reshape(ET0,[ET0.shape[0]/12,12,ET0.shape[1],ET0.shape[2]]); #full_ETs = np.zeros([12,180,360]); 
full_ETs = np.reshape(ETs,[ETs.shape[0]/12,12,ETs.shape[1],ETs.shape[2]]); #full_ETs = np.zeros([12,180,360]); 

#composite across EN or LN years
ENs1=w_s1[EN_index,...];ENs1=np.append(ENs1,EN_blank,axis=0)
ENs1=np.reshape(ENs1,[ENs1.shape[0]/42,42,ENs1.shape[1],ENs1.shape[2]])
ENs2=w_s2[EN_index,...];ENs2=np.append(ENs2,EN_blank,axis=0)
ENs2=np.reshape(ENs2,[ENs2.shape[0]/42,42,ENs2.shape[1],ENs2.shape[2]])
ENs3=w_s3[EN_index,...];ENs3=np.append(ENs3,EN_blank,axis=0)
ENs3=np.reshape(ENs3,[ENs3.shape[0]/42,42,ENs3.shape[1],ENs3.shape[2]])
ENcer=w_cer[EN_index,...];ENcer=np.append(ENcer,EN_blank,axis=0)
ENcer=np.reshape(ENcer,[ENcer.shape[0]/42,42,ENcer.shape[1],ENcer.shape[2]])
EN_ETcer=Ecer[EN_index,...];EN_ETcer=np.append(EN_ETcer,EN_blank,axis=0)
EN_ETcer=np.reshape(EN_ETcer,[EN_ETcer.shape[0]/42,42,EN_ETcer.shape[1],EN_ETcer.shape[2]])
EN_ETs=ETs[EN_index,...];EN_ETs=np.append(EN_ETs,EN_blank,axis=0)
EN_ETs=np.reshape(EN_ETs,[EN_ETs.shape[0]/42,42,EN_ETs.shape[1],EN_ETs.shape[2]])
EN_Pcer=p[EN_index,...];EN_Pcer=np.append(EN_Pcer,EN_blank,axis=0)
EN_Pcer=np.reshape(EN_Pcer,[EN_Pcer.shape[0]/42,42,EN_Pcer.shape[1],EN_Pcer.shape[2]])


LNs1=w_s1[LN_index,...];LNs1=np.append(LNs1,LN_blank,axis=0)
LNs1=np.reshape(LNs1,[LNs1.shape[0]/42,42,LNs1.shape[1],LNs1.shape[2]])
LNs2=w_s2[LN_index,...];LNs2=np.append(LNs2,LN_blank,axis=0)
LNs2=np.reshape(LNs2,[LNs2.shape[0]/42,42,LNs2.shape[1],LNs2.shape[2]])
LNs3=w_s3[LN_index,...];LNs3=np.append(LNs3,LN_blank,axis=0)
LNs3=np.reshape(LNs3,[LNs3.shape[0]/42,42,LNs3.shape[1],LNs3.shape[2]])
LNcer=w_cer[LN_index,...];LNcer=np.append(LNcer,LN_blank,axis=0)
LNcer=np.reshape(LNcer,[LNcer.shape[0]/42,42,LNcer.shape[1],LNcer.shape[2]])
LN_ETcer=Ecer[LN_index,...];LN_ETcer=np.append(LN_ETcer,LN_blank,axis=0)
LN_ETcer=np.reshape(LN_ETcer,[LN_ETcer.shape[0]/42,42,LN_ETcer.shape[1],LN_ETcer.shape[2]])
LN_ETs=ETs[LN_index,...];LN_ETs=np.append(LN_ETs,LN_blank,axis=0)
LN_ETs=np.reshape(LN_ETs,[LN_ETs.shape[0]/42,42,LN_ETs.shape[1],LN_ETs.shape[2]])
LN_Pcer=p[LN_index,...];LN_Pcer=np.append(LN_Pcer,LN_blank,axis=0)
LN_Pcer=np.reshape(LN_Pcer,[LN_Pcer.shape[0]/42,42,LN_Pcer.shape[1],LN_Pcer.shape[2]])

#general water stress timeseries.
ws1TS = np.nanmean(np.nanmean(w_s1[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)
ws2TS = np.nanmean(np.nanmean(w_s2[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)
ws3TS = np.nanmean(np.nanmean(w_s3[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)
wcerTS = np.nanmean(np.nanmean(w_cer[:,latMin:latMax,lonMin:lonMax],axis=2),axis=1)

#write the Time Series objects out
np.savetxt('/Volumes/Data_Archive/Results/soilWaterBalance/south1TS'+reg+anom+percent+'_lag'+str(lag)+str(years[0])+'_'+str(years[1])+'.txt',ws1TS)
np.savetxt('/Volumes/Data_Archive/Results/soilWaterBalance/south2TS'+reg+anom+percent+'_lag'+str(lag)+str(years[0])+'_'+str(years[1])+'.txt',ws2TS)
np.savetxt('/Volumes/Data_Archive/Results/soilWaterBalance/south3TS'+reg+anom+percent+'_lag'+str(lag)+str(years[0])+'_'+str(years[1])+'.txt',ws3TS)
np.savetxt('/Volumes/Data_Archive/Results/soilWaterBalance/cerradoTS'+reg+anom+percent+'_lag'+str(lag)+str(years[0])+'_'+str(years[1])+'.txt',wcerTS)


full_s1TS = np.nanmean(np.nanmean(full_s1[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
full_s2TS = np.nanmean(np.nanmean(full_s2[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
full_s3TS = np.nanmean(np.nanmean(full_s3[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
full_sTS = (full_s1TS+full_s2TS+full_s3TS)/3
full_cerTS = np.nanmean(np.nanmean(full_cer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
full_pTS = np.nanmean(np.nanmean(full_P[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
full_ETcerTS = np.nanmean(np.nanmean(full_ETcer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
full_ET0TS = np.nanmean(np.nanmean(full_ET0[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
full_ETsTS = np.nanmean(np.nanmean(full_ETs[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)


#full_s1TS = np.append(full_s1TS[:,7:],full_s1TS[:,:7],axis=1)
#full_s2TS = np.append(full_s2TS[:,7:],full_s2TS[:,:7],axis=1)
#full_s3TS = np.append(full_s3TS[:,7:],full_s3TS[:,:7],axis=1)
#full_sTS = np.append(full_sTS[:,7:],full_sTS[:,:7],axis=1)
#full_cerTS = np.append(full_cerTS[:,7:],full_cerTS[:,:7],axis=1)
#full_ETsTS = np.append(full_ETsTS[:,7:],full_ETsTS[:,:7],axis=1)
#full_ET0TS = np.append(full_ET0TS[:,7:],full_ET0TS[:,:7],axis=1)
#full_ETcerTS = np.append(full_ETcerTS[:,7:],full_ETcerTS[:,:7],axis=1)
#full_pTS = np.append(full_pTS[:,7:],full_pTS[:,:7],axis=1)

ENs1TS = np.nanmean(np.nanmean(ENs1[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
ENs2TS = np.nanmean(np.nanmean(ENs2[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
ENs3TS = np.nanmean(np.nanmean(ENs3[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
ENcerTS = np.nanmean(np.nanmean(ENcer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
EN_ETcerTS = np.nanmean(np.nanmean(EN_ETcer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
EN_ETsTS = np.nanmean(np.nanmean(EN_ETs[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
EN_PcerTS = np.nanmean(np.nanmean(EN_Pcer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)

LNs1TS = np.nanmean(np.nanmean(LNs1[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
LNs2TS = np.nanmean(np.nanmean(LNs2[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
LNs3TS = np.nanmean(np.nanmean(LNs3[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
LNcerTS = np.nanmean(np.nanmean(LNcer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
LN_ETcerTS = np.nanmean(np.nanmean(LN_ETcer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
LN_ETsTS = np.nanmean(np.nanmean(LN_ETs[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)
LN_PcerTS = np.nanmean(np.nanmean(LN_Pcer[:,:,latMin:latMax,lonMin:lonMax],axis=2),axis=2)


ENsTS = np.nanmean([ENs1TS[np.newaxis,...],ENs2TS[np.newaxis,...],ENs3TS[np.newaxis,...]],0)
LNsTS = np.nanmean([LNs1TS[np.newaxis,...],LNs2TS[np.newaxis,...],LNs3TS[np.newaxis,...]],0)
LNsTS = np.squeeze(LNsTS);ENsTS = np.squeeze(ENsTS)


full_s1Avg = np.nanmean(full_s1TS,axis=0);full_s2Avg = np.nanmean(full_s2TS,axis=0);
full_s3Avg = np.nanmean(full_s3TS,axis=0);full_cerAvg = np.nanmean(full_cerTS,axis=0);
full_ETcerAvg = np.nanmean(full_ETcerTS,axis=0);full_ET0Avg = np.nanmean(full_ET0TS,axis=0);
full_pAvg = np.nanmean(full_pTS,axis=0);full_sAvg = np.nanmean(full_sTS,axis=0)
full_s1Min = np.nanmin(full_s1TS,axis=0);full_s2Min = np.nanmin(full_s2TS,axis=0);
full_s3Min = np.nanmin(full_s3TS,axis=0);full_cerMin = np.nanmin(full_cerTS,axis=0);
full_ETcerMin = np.nanmin(full_ETcerTS,axis=0);full_ET0Min = np.nanmin(full_ET0TS,axis=0);
full_pMin = np.nanmin(full_pTS,axis=0);full_sMin = np.nanmin(full_sTS,axis=0)
full_s1Max = np.nanmax(full_s1TS,axis=0);full_s2Max = np.nanmax(full_s2TS,axis=0);
full_s3Max = np.nanmax(full_s3TS,axis=0);full_cerMax = np.nanmax(full_cerTS,axis=0);
full_ETcerMax = np.nanmax(full_ETcerTS,axis=0);full_ET0Max = np.nanmax(full_ET0TS,axis=0);
full_pMax = np.nanmax(full_pTS,axis=0);full_sMax = np.nanmax(full_sTS,axis=0)

full_ETsAvg = np.nanmean(full_ETsTS,axis=0)
full_ETsMin = np.nanmin(full_ETsTS,axis=0)
full_ETsMax = np.nanmax(full_ETsTS,axis=0)

ENsTSavg = np.nanmean(ENsTS,axis=0);LNsTSavg = np.nanmean(LNsTS,axis=0)
ENs1TSavg = np.nanmean(ENs1TS,axis=0);LNs1TSavg = np.nanmean(LNs1TS,axis=0)
ENs2TSavg = np.nanmean(ENs2TS,axis=0);LNs2TSavg = np.nanmean(LNs2TS,axis=0)
ENs3TSavg = np.nanmean(ENs3TS,axis=0);LNs3TSavg = np.nanmean(LNs3TS,axis=0)
ENcerTSavg = np.nanmean(ENcerTS,axis=0);LNcerTSavg = np.nanmean(LNcerTS,axis=0)
EN_ETcerTSavg = np.nanmean(EN_ETcerTS,axis=0);LN_ETcerTSavg = np.nanmean(LN_ETcerTS,axis=0)
EN_PcerTSavg = np.nanmean(EN_PcerTS,axis=0);LN_PcerTSavg = np.nanmean(LN_PcerTS,axis=0)

ENsTSmin = np.nanmin(ENsTS,axis=0);LNsTSmin = np.nanmin(LNsTS,axis=0)
ENs1TSmin = np.nanmin(ENs1TS,axis=0);LNs1TSmin = np.nanmin(LNs1TS,axis=0)
ENs2TSmin = np.nanmin(ENs2TS,axis=0);LNs2TSmin = np.nanmin(LNs2TS,axis=0)
ENs3TSmin = np.nanmin(ENs3TS,axis=0);LNs3TSmin = np.nanmin(LNs3TS,axis=0)
ENcerTSmin = np.nanmin(ENcerTS,axis=0);LNcerTSmin = np.nanmin(LNcerTS,axis=0)
EN_ETcerTSmin = np.nanmin(EN_ETcerTS,axis=0);LN_ETcerTSmin = np.nanmin(LN_ETcerTS,axis=0)
EN_PcerTSmin = np.nanmin(EN_PcerTS,axis=0);LN_PcerTSmin = np.nanmin(LN_PcerTS,axis=0)

ENsTSmax = np.nanmax(ENsTS,axis=0);LNsTSmax = np.nanmax(LNsTS,axis=0)
ENs1TSmax = np.nanmax(ENs1TS,axis=0);LNs1TSmax = np.nanmax(LNs1TS,axis=0)
ENs2TSmax = np.nanmax(ENs2TS,axis=0);LNs2TSmax = np.nanmax(LNs2TS,axis=0)
ENs3TSmax = np.nanmax(ENs3TS,axis=0);LNs3TSmax = np.nanmax(LNs3TS,axis=0)
ENcerTSmax = np.nanmax(ENcerTS,axis=0);LNcerTSmax = np.nanmax(LNcerTS,axis=0)
EN_ETcerTSmax = np.nanmax(EN_ETcerTS,axis=0);LN_ETcerTSmax = np.nanmax(LN_ETcerTS,axis=0)
EN_PcerTSmax = np.nanmax(EN_PcerTS,axis=0);LN_PcerTSmax = np.nanmax(LN_PcerTS,axis=0)

#monNums = [5,9,13,17,21,25,29,33,37,41]
#monNams = ['','Oct','Feb','Jun','Oct','Feb','Jun','Oct','Feb','Jun']
monNums = [5,8,13,17,20,25,29,32,37,41]
monNams = ['','Sep','Feb','Jun','Sep','Feb','Jun','Sep','Feb','Jun']

if plot is 'SWC':
    fig = plt.figure();
    ax1=plt.subplot(321);ax2=plt.subplot(323);ax3=plt.subplot(325)
    ax4=plt.subplot(322);ax5=plt.subplot(324);ax6=plt.subplot(326)
    ax1.set_title(u'El Niño\n                  EN -1                     EN 0                       EN +1');
    ax4.set_title(u'La Niña\n                  LN -1                     LN 0                       LN +1');
    ax1.plot(range(42),EN_PcerTSavg[:42],'k',lw=4); ax1.plot(EN_PcerTS.T,ls='--',color='k',alpha=0.25);
    ax1.set_ylim([40,240]);ax1.set_xlim([0,41]);ax1.set_ylabel('Precipitation')
    ax2.plot(range(42),np.nanmean(ENcerTS[:,:42].T,axis=1),'k',lw=3); ax2.plot(ENcerTS.T,color='k',ls='--',alpha=0.3);
    ax3.plot(range(42),np.nanmean(ENsTS[:,:42].T,axis=1),'k',lw=3); ax3.plot(ENsTS.T,'--k',alpha=0.5);
    ax2.set_ylim([.6,1]);ax3.set_ylim([.6,1]);ax2.set_xlim([0,41]);ax3.set_xlim([0,41]);
    ax2.set_ylabel('Safrinha \n% Soil Water Content');ax3.set_ylabel('Early Season \n% Soil Water Content')
    
    #ax1.axvspan(13,17,facecolor='cornflowerblue', alpha=0.20);ax1.axvspan(25,29,facecolor='cornflowerblue', alpha=0.20);ax1.axvspan(37,41,facecolor='cornflowerblue', alpha=0.20);
    #ax1.axvspan(9,16,facecolor='k', alpha=0.20);ax1.axvspan(21,28,facecolor='k', alpha=0.20);ax1.axvspan(33,40,facecolor='k', alpha=0.20);
    #ax1.axvspan(8,13,facecolor='b', alpha=0.20);ax1.axvspan(20,25,facecolor='b', alpha=0.20);ax1.axvspan(32,37,facecolor='b', alpha=0.20);
    ax2.axvspan(13,17,facecolor='cornflowerblue', alpha=0.20);ax2.axvspan(25,29,facecolor='cornflowerblue', alpha=0.20);ax2.axvspan(37,41,facecolor='cornflowerblue', alpha=0.20);
    ax3.axvspan(9,16,facecolor='darkgoldenrod', alpha=0.40);ax3.axvspan(21,28,facecolor='darkgoldenrod', alpha=0.40);ax3.axvspan(33,40,facecolor='darkgoldenrod', alpha=0.40);
    ax2.axvspan(8,13,facecolor='b', alpha=0.20);ax2.axvspan(20,25,facecolor='b', alpha=0.20);ax2.axvspan(32,37,facecolor='b', alpha=0.20);
    ax1.set_xticks(monNums);ax1.set_xticklabels(monNams);
    ax2.set_xticks(monNums);ax2.set_xticklabels(monNams);
    ax3.set_xticks(monNums);ax3.set_xticklabels(monNams);
    
    ax4.plot(range(42),LN_PcerTSavg[:42],'k',lw=4); ax4.plot(LN_PcerTS.T,ls='--',color='k',alpha=0.25);
    ax4.set_ylim([40,240]);ax4.set_xlim([0,41]);#ax4.set_ylabel('Precipitation')
    ax5.plot(range(42),np.nanmean(LNcerTS[:,:42].T,axis=1),'k',lw=3); ax5.plot(LNcerTS[:,:42].T,ls='--',color='k',alpha=0.4);
    ax6.plot(range(42),LNsTSavg[:42].T,'k',lw=3); ax6.plot(LNsTS.T,'--k',alpha=0.5);
    ax5.set_ylim([.6,1]);ax5.set_xlim([0,41]);ax6.set_ylim([.6,1]);ax6.set_xlim([0,41]);
    
    #ax4.axvspan(13,17,facecolor='cornflowerblue', alpha=0.20);ax4.axvspan(25,29,facecolor='cornflowerblue', alpha=0.20);ax4.axvspan(37,41,facecolor='cornflowerblue', alpha=0.20);
    #ax4.axvspan(8,13,facecolor='b', alpha=0.20);ax4.axvspan(20,25,facecolor='b', alpha=0.20);ax4.axvspan(32,37,facecolor='b', alpha=0.20);
    #ax4.axvspan(9,16,facecolor='k', alpha=0.20);ax4.axvspan(21,28,facecolor='k', alpha=0.20);ax4.axvspan(33,40,facecolor='k', alpha=0.20);
    ax5.axvspan(8,13,facecolor='b', alpha=0.20);ax5.axvspan(20,25,facecolor='b', alpha=0.20);ax5.axvspan(32,37,facecolor='b', alpha=0.20);
    ax5.axvspan(13,17,facecolor='cornflowerblue', alpha=0.20);ax5.axvspan(25,29,facecolor='cornflowerblue', alpha=0.20);ax5.axvspan(37,41,facecolor='cornflowerblue', alpha=0.20);
    ax6.axvspan(9,16,facecolor='darkgoldenrod', alpha=0.40);ax6.axvspan(21,28,facecolor='darkgoldenrod', alpha=0.40);ax6.axvspan(33,40,facecolor='darkgoldenrod', alpha=0.40);
    ax4.set_xticks(monNums);ax4.set_xticklabels(monNams);
    ax5.set_xticks(monNums);ax5.set_xticklabels(monNams)
    ax6.set_xticks(monNums);ax6.set_xticklabels(monNams)
    
    fig.set_size_inches(16, 10)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/SWC_ENSO/SWCenso_'+reg+'.png')
    plt.close()


if plot is 'P_tot':
    ############# Plot of LN cycle for moisture, E and P ################
    fig = plt.figure();
    ax1=plt.subplot(321);ax2=plt.subplot(323);ax3=plt.subplot(325)
    ax4=plt.subplot(322);ax5=plt.subplot(324);ax6=plt.subplot(326)
    ax1.set_title(u'El Niño\n');ax4.set_title(u'La Niña\n')
    ax1.plot(range(42),EN_PcerTSavg[:42],'b',lw=4); ax1.plot(EN_PcerTS.T,'--b',alpha=0.25);
    ax1.set_ylim([30,250]);ax1.set_xlim([0,41]);ax1.set_ylabel('Precipitation')
    ax2.plot(range(42),np.nanmean(EN_ETcerTS[:,:42].T,axis=1),'g',lw=3); ax2.plot(EN_ETcerTS[:,:42].T,'--g',alpha=0.25);
    ax2.plot(range(42),np.nanmean(EN_ETsTS[:,:42].T,axis=1),'k',lw=3); ax2.plot(EN_ETsTS[:,:42].T,'--k',alpha=0.25);
    ax2.set_ylim([0,150]);ax2.set_xlim([0,41]);ax2.set_ylabel('Evapotranspiration')
    ax3.plot(range(42),np.nanmean(ENcerTS[:,:42].T,axis=1),'g',lw=3); ax3.plot(ENcerTS.T,'--g',alpha=0.25);
    ax3.plot(range(42),np.nanmean(ENsTS[:,:42].T,axis=1),'k',lw=3); ax3.plot(ENsTS.T,'--k',alpha=0.25);
    ax3.set_ylim([.5,1]);ax3.set_xlim([0,41]);ax3.set_ylabel('% Soil Water Content')
#    ax1.axvspan(13,17,facecolor='g', alpha=0.20);ax1.axvspan(25,29,facecolor='g', alpha=0.20);ax1.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax1.axvspan(10,16,facecolor='k', alpha=0.20);ax1.axvspan(22,28,facecolor='k', alpha=0.20);ax1.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax1.axvspan(9,13,facecolor='y', alpha=0.20);ax1.axvspan(21,25,facecolor='y', alpha=0.20);ax1.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax2.axvspan(13,17,facecolor='g', alpha=0.20);ax2.axvspan(25,29,facecolor='g', alpha=0.20);ax2.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax2.axvspan(10,16,facecolor='k', alpha=0.20);ax2.axvspan(22,28,facecolor='k', alpha=0.20);ax2.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax2.axvspan(9,13,facecolor='y', alpha=0.20);ax2.axvspan(21,25,facecolor='y', alpha=0.20);ax2.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax3.axvspan(13,17,facecolor='g', alpha=0.20);ax3.axvspan(25,29,facecolor='g', alpha=0.20);ax3.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax3.axvspan(10,16,facecolor='k', alpha=0.20);ax3.axvspan(22,28,facecolor='k', alpha=0.20);ax3.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax3.axvspan(9,13,facecolor='y', alpha=0.20);ax3.axvspan(21,25,facecolor='y', alpha=0.20);ax3.axvspan(33,37,facecolor='y', alpha=0.20);
    
    ax4.plot(range(42),LN_PcerTSavg[:42],'b',lw=4); ax4.plot(LN_PcerTS.T,'--b',alpha=0.25);
    ax4.set_ylim([30,250]);ax4.set_xlim([0,41]);#ax4.set_ylabel('Precipitation')
    ax5.plot(range(42),np.nanmean(LN_ETcerTS[:,:42].T,axis=1),'g',lw=3); ax5.plot(LN_ETcerTS[:,:42].T,'--g',alpha=0.25);
    ax5.plot(range(42),np.nanmean(LN_ETsTS[:,:42].T,axis=1),'k',lw=3); ax5.plot(LN_ETsTS[:,:42].T,'--k',alpha=0.25);
    ax5.set_ylim([0,150]);ax5.set_xlim([0,41]);#ax5.set_ylabel('Evapotranspiration')
    ax6.plot(range(42),np.nanmean(LNcerTS[:,:42].T,axis=1),'g',lw=3); ax6.plot(LNcerTS[:,:42].T,'--g',alpha=0.25);
    ax6.plot(range(42),LNsTSavg[:42].T,'k',lw=3); ax6.plot(LNsTS.T,'--k',alpha=0.25);
    ax6.set_ylim([.5,1]);ax6.set_xlim([0,41]);#ax6.set_ylabel('% Soil Water Content')
#    ax4.axvspan(13,17,facecolor='g', alpha=0.20);ax4.axvspan(25,29,facecolor='g', alpha=0.20);ax4.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax4.axvspan(9,13,facecolor='y', alpha=0.20);ax4.axvspan(21,25,facecolor='y', alpha=0.20);ax4.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax4.axvspan(10,16,facecolor='k', alpha=0.20);ax4.axvspan(22,28,facecolor='k', alpha=0.20);ax4.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax5.axvspan(9,13,facecolor='y', alpha=0.20);ax5.axvspan(21,25,facecolor='y', alpha=0.20);ax5.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax5.axvspan(13,17,facecolor='g', alpha=0.20);ax5.axvspan(25,29,facecolor='g', alpha=0.20);ax5.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax5.axvspan(10,16,facecolor='k', alpha=0.20);ax5.axvspan(22,28,facecolor='k', alpha=0.20);ax5.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax6.axvspan(9,13,facecolor='y', alpha=0.20);ax6.axvspan(21,25,facecolor='y', alpha=0.20);ax6.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax6.axvspan(13,17,facecolor='g', alpha=0.20);ax6.axvspan(25,29,facecolor='g', alpha=0.20);ax6.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax6.axvspan(10,16,facecolor='k', alpha=0.20);ax6.axvspan(22,28,facecolor='k', alpha=0.20);ax6.axvspan(34,40,facecolor='k', alpha=0.20);
    

if plot is 'P_anom':
    ############ Plot Anoms ################
    fig = plt.figure();fig.suptitle(u'La Niña')
    ax1=plt.subplot(311);ax2=plt.subplot(312);ax3=plt.subplot(313)
    ax1.plot(range(42),(LN_PcerTSavg[:42]-np.tile(full_pAvg,4)[:42]),'b',lw=4); #ax1.plot(LN_PcerTS.T-np.repeat(np.tile(full_pAvg,4)[:42,np.newaxis],10,1),'--b',alpha=0.25);
    ax1.set_ylim([-35,35]);ax1.set_xlim([0,41]);ax1.set_ylabel('Precipitation')
    ax2.plot(range(42),(np.nanmean(LN_ETcerTS[:,:42],axis=0)-np.tile(full_ETcerAvg,4)[:42]).T,'g',lw=3);# ax2.plot(LN_ETcerTS[:,:42].T,'--g',alpha=0.25);
    ax2.plot(range(42),(np.nanmean(LN_ETsTS[:,:42],axis=0)-np.tile(full_ETsAvg,4)[:42]).T,'k',lw=3);# ax2.plot(LN_ETsTS[:,:42].T,'--k',alpha=0.25);
    ax2.set_ylim([-10,10]);ax2.set_xlim([0,41]);ax2.set_ylabel('Evapotranspiration')
    ax3.plot(range(42),(LNcerTSavg[:42]-np.tile(full_cerAvg,4)[:42]).T,'g',lw=3);# ax3.plot(LNcerTS[:,:42].T,'--g',alpha=0.25);
    ax3.plot(range(42),(LNsTSavg[:42]-np.tile(full_sAvg,4)[:42]).T,'k',lw=3);# ax3.plot(LNsTS.T,'--k',alpha=0.25);
    ax3.set_ylim([-.10,.10]);ax3.set_xlim([0,41]);ax3.set_ylabel('% Soil Water Content')
    ax1.hlines(0,0,42);ax2.hlines(0,0,42);ax3.hlines(0,0,42)
#    ax1.axvspan(13,17,facecolor='g', alpha=0.20);ax1.axvspan(25,29,facecolor='g', alpha=0.20);ax1.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax1.axvspan(9,13,facecolor='y', alpha=0.20);ax1.axvspan(21,25,facecolor='y', alpha=0.20);ax1.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax1.axvspan(10,16,facecolor='k', alpha=0.20);ax1.axvspan(22,28,facecolor='k', alpha=0.20);ax1.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax2.axvspan(9,13,facecolor='y', alpha=0.20);ax2.axvspan(21,25,facecolor='y', alpha=0.20);ax2.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax2.axvspan(13,17,facecolor='g', alpha=0.20);ax2.axvspan(25,29,facecolor='g', alpha=0.20);ax2.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax2.axvspan(10,16,facecolor='k', alpha=0.20);ax2.axvspan(22,28,facecolor='k', alpha=0.20);ax2.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax3.axvspan(9,13,facecolor='y', alpha=0.20);ax3.axvspan(21,25,facecolor='y', alpha=0.20);ax3.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax3.axvspan(13,17,facecolor='g', alpha=0.20);ax3.axvspan(25,29,facecolor='g', alpha=0.20);ax3.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax3.axvspan(10,16,facecolor='k', alpha=0.20);ax3.axvspan(22,28,facecolor='k', alpha=0.20);ax3.axvspan(34,40,facecolor='k', alpha=0.20);
  
    
    fig = plt.figure();fig.suptitle(u'El Niño')
    ax1=plt.subplot(311);ax2=plt.subplot(312);ax3=plt.subplot(313)
    ax1.plot(range(42),(EN_PcerTSavg[:42]-np.tile(full_pAvg,4)[:42]),'b',lw=4); #ax1.plot(EN_PcerTS.T-np.repeat(np.tile(full_pAvg,4)[:42,np.newaxis],10,1),'--b',alpha=0.25);
    ax1.set_ylim([-35,35]);ax1.set_xlim([0,41]);ax1.set_ylabel('Precipitation')
    ax2.plot(range(42),(np.nanmean(EN_ETcerTS[:,:42],axis=0)-np.tile(full_ETcerAvg,4)[:42]).T,'g',lw=3);# ax2.plot(EN_ETcerTS[:,:42].T,'--g',alpha=0.25);
    ax2.plot(range(42),(np.nanmean(EN_ETsTS[:,:42],axis=0)-np.tile(full_ETsAvg,4)[:42]).T,'k',lw=3);# ax2.plot(EN_ETsTS[:,:42].T,'--k',alpha=0.25);
    ax2.set_ylim([-20,20]);ax2.set_xlim([0,41]);ax2.set_ylabel('Evapotranspiration')
    ax3.plot(range(42),(ENcerTSavg[:42]-np.tile(full_cerAvg,4)[:42]).T,'g',lw=3);# ax3.plot(ENcerTS[:,:42].T,'--g',alpha=0.25);
    ax3.plot(range(42),(ENsTSavg[:42]-np.tile(full_sAvg,4)[:42]).T,'k',lw=3);# ax3.plot(ENsTS.T,'--k',alpha=0.25);
    ax3.set_ylim([-.10,.10]);ax3.set_xlim([0,41]);ax3.set_ylabel('% Soil Water Content')
    ax1.hlines(0,0,42);ax2.hlines(0,0,42);ax3.hlines(0,0,42)
#    ax1.axvspan(13,17,facecolor='g', alpha=0.20);ax1.axvspan(25,29,facecolor='g', alpha=0.20);ax1.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax1.axvspan(9,13,facecolor='y', alpha=0.20);ax1.axvspan(21,25,facecolor='y', alpha=0.20);ax1.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax1.axvspan(10,16,facecolor='k', alpha=0.20);ax1.axvspan(22,28,facecolor='k', alpha=0.20);ax1.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax2.axvspan(9,13,facecolor='y', alpha=0.20);ax2.axvspan(21,25,facecolor='y', alpha=0.20);ax2.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax2.axvspan(13,17,facecolor='g', alpha=0.20);ax2.axvspan(25,29,facecolor='g', alpha=0.20);ax2.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax2.axvspan(10,16,facecolor='k', alpha=0.20);ax2.axvspan(22,28,facecolor='k', alpha=0.20);ax2.axvspan(34,40,facecolor='k', alpha=0.20);
#    ax3.axvspan(9,13,facecolor='y', alpha=0.20);ax3.axvspan(21,25,facecolor='y', alpha=0.20);ax3.axvspan(33,37,facecolor='y', alpha=0.20);
#    ax3.axvspan(13,17,facecolor='g', alpha=0.20);ax3.axvspan(25,29,facecolor='g', alpha=0.20);ax3.axvspan(37,41,facecolor='g', alpha=0.20);
#    ax3.axvspan(10,16,facecolor='k', alpha=0.20);ax3.axvspan(22,28,facecolor='k', alpha=0.20);ax3.axvspan(34,40,facecolor='k', alpha=0.20);
    
    

#plt.plot(range(42),ENcerTSavg,'k');plt.fill_between(range(42),ENcerTSmin,ENcerTSmax,alpha=0.25,color='k')
#plt.plot(range(42),EN_ETcerTSavg,'g');plt.fill_between(range(42),EN_ETcerTSmin,EN_ETcerTSmax,alpha=0.25,color='g')
#plt.plot(range(42),EN_PcerTSavg,'b');plt.fill_between(range(42),EN_PcerTSmin,EN_PcerTSmax,alpha=0.25,color='b')


############# Time Series plots ##############
#La Niña soil water balance plots
#plt.plot(range(42),LNcerTSavg,'k');plt.fill_between(range(42),LNcerTSmin,LNcerTSmax,alpha=0.25,color='k')
#plt.plot(range(42),LNsTSavg,'b');plt.fill_between(range(42),LNsTSmin,LNsTSmax,alpha=0.25,color='b')
#plt.plot(range(42),LNs1TSavg,'b');plt.fill_between(range(42),LNs1TSmin,LNs1TSmax,alpha=0.25,color='b')
#plt.plot(range(42),LNs2TSavg,'b');plt.fill_between(range(42),LNs2TSmin,LNs2TSmax,alpha=0.25,color='b')
#plt.plot(range(42),LNs3TSavg,'b');plt.fill_between(range(42),LNs3TSmin,LNs3TSmax,alpha=0.25,color='b')

#El Niño soil water balance plots
#plt.plot(range(42),ENcerTSavg,'k');plt.fill_between(range(42),ENcerTSmin,ENcerTSmax,alpha=0.25,color='k')
#plt.plot(range(42),ENsTSavg,'b');plt.fill_between(range(42),ENsTSmin,ENsTSmax,alpha=0.25,color='b')
#plt.plot(range(42),ENs1TSavg,'b');plt.fill_between(range(42),ENs1TSmin,ENs1TSmax,alpha=0.25,color='b')
#plt.plot(range(42),ENs2TSavg,'b');plt.fill_between(range(42),ENs2TSmin,ENs2TSmax,alpha=0.25,color='b')
#plt.plot(range(42),ENs3TSavg,'b');plt.fill_between(range(42),ENs3TSmin,ENs3TSmax,alpha=0.25,color='b')


