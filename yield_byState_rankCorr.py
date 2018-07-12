# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 13:57:51 2015

@author: weston

Sep 29 2015 - updated to include CONAB results in Brazil for soybean
Oct 21 2015 - updated soybean anomalies to include first difference and % change methods
"""

import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy import signal
import netCDF4
from mpl_toolkits.basemap import Basemap
from scipy.stats.stats import spearmanr
start = time.clock()

countries = ['BR']
proj = 'AmENSO_yields' #AusENSO or AmENSO_yields or ChnENSO
printState = 'N'#Y/N
corrThresh = .1
numSurrogates = 1000
stORreg = 'state' #state or region
BRmaize = 'maize' 
#NOTE: CHANGE TO SPRING OR WINTER WHEAT FURTHER DOWN IN THE "US" SECTION

smooth = 4
def moving_average(a, n=smooth*2+1) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
  
#create surrogate timeseries
#follows from Ebusuzaki (1997) and Schrieber and Shmitz (2000)
def surrogate(ai):
    ak = np.fft.fft(ai) #DFFT
    k = np.size(ai)
    rk = [None]*k
    rk[0] = ak[0]*np.sin(ak[0]) #amp zero frequency (with sign)
    if k%2==0: #if even
        for i in range(1,k//2):
            thetak = np.pi*2.*np.random.uniform(0,1)   #generate random phase      
            rk[i] = np.abs(ak[i])*np.exp(thetak*1j)    #apply random phase to positive freq
            rk[-(i)] = np.abs(ak[-(i+1)])*np.exp(thetak*-1j) #apply random phase to negative freq
        thetak = np.pi*2.*np.random.uniform(0,1) #generate random phase    
        rk[k//2]=np.sqrt(2.)*np.abs(ak[k//2])*np.cos(thetak) #amp the nyquist freq
    else: #if odd
        for i in range(1,(k+1)//2):
            thetak = np.pi*2.*np.random.uniform(0,1)      #generate random phase       
            rk[i] = np.abs(ak[i])*np.exp(thetak*1j) #apply random phase to positive freq
            rk[-(i)] = np.abs(ak[i])*np.exp(thetak*-1j) #apply random phase to negative freq
    rj = np.fft.ifft(rk) #inverse fft
    return rj.real

for country in countries:

    if country is 'BR':
        lag = ['N','N','Y']
        dataFile = 'BR_CONAB'#'BR_IBGE', 'BR_CONAB'
        stYr = [1978,1978,1978]#, [1991,1990,1991]#
        endYr = [2013,2013,2013]
        mon= [10,10,10]#[10,6,6] #[wheat, maize, soy], NOTE PYTHON INDEXING
        ensoMon = [10,10,10]#[10,6,6] #this always starts at 1950, so have to use 0 for previous in US and lag as 'N'
        crop = 'maize'
        if (BRmaize=='maize1')|(BRmaize=='maize'):
            breakPt = 12 #number of years in to put the breakpoint for detrending
        elif BRmaize=='maize2':
            stYr = [1991,1990,1991]#
            breakPt = 0
    elif country is 'AR':
        lag = ['N','N','N'] #BR has to be lagged for soy to match AR. Both reflect planting year then
        dataFile = 'AR_SIIA'
        stYr = [1969,1969,1970]
        endYr = [2012,2012,2013]
        mon= [10,10,10]#[10,6,6] #[wheat, maize, soy], NOTE PYTHON INDEXING
        ensoMon = [10,10,10]#[10,6,6] #this always starts at 1950, so have to use 0 for previous in US and lag as 'N'
        breakPt = 0#25 #number of years in to put the breakpoint for detrending
        crop = 'all'
    elif country is 'UY':
        lag = ['N','N','N']
        dataFile = 'UY'
        stYr = [1985,1985,1985]
        endYr = [1997,2009,1997]
        mon= [9,11,10]#[10,6,6] #[wheat, maize, soy], NOTE PYTHON INDEXING
        ensoMon = [9,11,1]#[10,6,6] #this always starts at 1950, so have to use 0 for previous in US and lag as 'N'
        breakPt = 0 #number of years in to put the breakpoint for detrending
        crop = 'all'
    elif country is 'US':
        lag = ['Y','N','N']
        dataFile = 'USDA_NASS'
        stYr = [1950,1950,1950]
        endYr = [2012,2012,2012]
        mon= [10,10,10] #[4,6,6]#[wheat, maize, soy], NOTE PYTHON INDEXING
        ensoMon = [10,10,10] #[4,6,6]#this always starts at 1950, so have to use 0 for previous in US and lag as 'N'
        breakPt = 0 #number of years in to put the breakpoint for detrending   
        crop = 'all'
        USw = 'winterwheat'
    elif country is 'MX':
        lag = ['N','N','N']
        dataFile = 'MX_SAGARPA'
        stYr = [1980,1980,1980]
        endYr = [2012,2012,2012]
        mon= [10,10,10]#[10,6,6] #[wheat, maize, soy], NOTE PYTHON INDEXING
        ensoMon = [10,10,10]#[10,6,6] #this always starts at 1950, so have to use 0 for previous in US and lag as 'N'
        breakPt = 0 #number of years in to put the breakpoint for detrending   
        crop = 'maize'
    elif country is 'CA':
        lag = ['N','N','N']
        dataFile = 'CA_CANSIM'
        stYr = [1950,1950,1950]
        endYr = [2012,2012,2012]
        mon= [6,6,6]#[10,6,6] #[wheat, maize, soy], NOTE PYTHON INDEXING
        ensoMon = [6,6,6]#[10,6,6] #this always starts at 1950, so have to use 0 for previous in US and lag as 'N'
        breakPt = 0 #number of years in to put the breakpoint for detrending   
        crop = 'wheat'
    elif country is 'AUS':
        lag = ['N','N','N']
        dataFile = 'AUS_ABS'
        stYr = [1950,1950,1950]
        endYr = [2012,2012,2012]
        mon= [8,6,6]# NOTE PYTHON INDEXING
        ensoMon = [8,6,6]#
        breakPt = 30 #number of years in to put the breakpoint for detrending   
        crop = 'wheat'
    elif country is 'CHN':
        lag = ['Y','N','N']
        dataFile = 'CHN'
        stYr = [1950,1950,1950]
        endYr = [2012,2012,2012]
        mon= [5,7,6]#NOTE PYTHON INDEXING
        ensoMon = [5,7,6]
        crop = 'maize'
        if(crop=='wheat'):
            mon= [5,7,6]#NOTE PYTHON INDEXING
            ensoMon = [5,7,6]
        elif(crop=='winterWheat'):
            crop='wheat'
            mon= [3,7,6]#NOTE PYTHON INDEXING
            ensoMon = [3,7,6]
        breakPt = 30 #number of years in to put the breakpoint for detrending   

    months=['January','February','March','April','May','June','July','August','September','October','Novemver','December']
    
    wLab = str(months[mon[0]])+' SST'
    mLab = str(months[mon[1]])+' SST'
    soLab = str(months[mon[2]])+' SST'
    
    yldSt = [0,0,0]
    yldEnd = [0,0,0]
    for iy in range(3):
        if lag[iy] is 'Y': yldSt[iy]=stYr[iy]+1; yldEnd[iy]=endYr[iy]
        elif lag[iy] is '-Y': 
            yldSt[iy]=stYr[iy]-1; yldEnd[iy]=endYr[iy]-2
        else: yldSt[iy]=stYr[iy]; yldEnd[iy]=endYr[iy]-1
    
    #Read in ENSO information 
    enso34tab = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/ENSO/NinoIndex34.csv', 
                         header=1)
    ensoAll = np.array(enso34tab.loc[str(stYr[0]):str(endYr[0]-1)])
    ensoAll = ensoAll[:,0:12]
    
    ensoTS = np.zeros([ensoAll.shape[0],3,numSurrogates])*np.nan #create the dataframe for the enso Surrogate time series
    
    #create the ensemble of ENSO random timeseries through phase randomization
    for im in range(3):
        for isu in range(numSurrogates):
            ensnoTS[:,im,isu]=surrogate(ensoAll[:,im])    
    
    if (crop != 'soy')&(crop != 'maize'):
        if country is 'US':
            USwh = USw
        else: USwh = 'wheat'
                    #---------------------------------- Wheat ------------------------------------#
        #read in the yield data
        CNTRYyld = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/AmENSO_yld/AmENSO_'+country+'_'+USwh+'_halfPercMaskGDD2all'+stORreg+'.csv', 
                         header=0)
        CNTRYyld.set_index([CNTRYyld.index.year,'state'],drop=True,inplace=True)           
        CNTRYyld = CNTRYyld['yldAnomGau'];CNTRYyld=CNTRYyld.sort_index()
#        CNTRYyld = CNTRYyld.sort()
        
        #read in the harvested area data
        CNTRYha = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/HarvestedArea/HistoricalData/'+dataFile+'/'+country+'_'+USwh+'_'+stORreg+'.csv',
                                      header=0)
        CNTRYha = CNTRYha.sort()
        
        #read in the total production data
        CNTRYp = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/production/HistoricalData/'+dataFile+'/'+country+'_'+USwh+'_'+stORreg+'.csv',
                                      header=0)
        CNTRYp = CNTRYp.sort(); #CNTRYp = np.array(CNTRYp[str(yldSt[0]):str(yldEnd[0])]['Value'])
#        if dataFile == 'BR_CONAB':
#            CNTRYp=CNTRYp*1000.
        #NaNs coded as 0 in the data
        CNTRYyld[CNTRYyld==0]=np.nan
        CNTRYha[CNTRYha==0]=np.nan
        
        #Read in the SST anomaly data
        sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
        #T is months since 1960-01-01, and I want to start at 1950 with one year lag
        sst=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-stYr[0]))&
                                    (sstData.variables['T'][:]<=12*(endYr[0]-1960)),:,
                                    (sstData.variables['Y'][:]<45)&
                                    (sstData.variables['Y'][:]>-45),:]
        #pull only november SST anomalies
        sst=sst[(np.array(range(endYr[0]-stYr[0]))*12)+mon[0],...]
        sstLats=sstData.variables['Y'][(sstData.variables['Y'][:]<45)&
                                    (sstData.variables['Y'][:]>-45)]; sstLons=sstData.variables['X'][:]
        sstLons,sstLats=np.meshgrid(sstLons,sstLats)
           
        REGp1 = []
        states1 = []
        statesR1 = []
        statesSig1 = []
        stateSST = np.zeros(np.shape(sstLons))
        
        for STnam in np.unique(CNTRYyld.index.get_level_values(1)):
            surrogateR1 = [] #clear previous surrogate correlations
            if np.array(CNTRYyld.loc[(yldSt[0]):(yldEnd[0]),STnam]).dtype == 'O':
                continue#this skips columns with missing values
            ST_temp = CNTRYyld.loc[(yldSt[0]):(yldEnd[0]),STnam]
            ensoAll_temp = ensoAll[np.isfinite(np.array(ST_temp)),ensoMon[0]]
            ensoTS_temp = ensoTS[np.isfinite(np.array(ST_temp)),0,:] #clip the surrogate TS to dates
#            ensoAll_temp =np.mean([ensoAll[np.isfinite(np.array(ST_temp)),ensoMon[0]]-1,
#                               ensoAll[np.isfinite(np.array(ST_temp)),ensoMon[0]],
#                               ensoAll[np.isfinite(np.array(ST_temp)),ensoMon[0]]+1],axis=0)
            sst_temp = sst[np.isfinite(np.array(ST_temp)),...]
            ST_temp = ST_temp.dropna()
            STyrs = ST_temp.index
            if len(ST_temp) < endYr[0]-stYr[0]-1:
                continue
            ST = ST_temp#np.array(signal.detrend(ST_temp,bp=breakPt)/(ST_temp-signal.detrend(ST_temp,bp=breakPt))) #%yield anom = (yld obs - yld exp)/yld exp
            for isu in range(numSurrogates): #correlate against surrogate TS to find significance
                surrogateR1.append(spearmanr(ensoTS_temp[:,isu],ST)[0])
                    #sort and find non parametric alphas
            surrogateR1 = np.sort(surrogateR1)
            LowBound = surrogateR1[numSurrogates*corrThresh//2]
            HighBound = surrogateR1[-numSurrogates*corrThresh//2]
            STr = spearmanr(ST,ensoAll_temp)[0]
            sig = spearmanr(ST,ensoAll_temp)[1]
            if printState is 'Y':
                for iy in range(sstLats.shape[0]):
                    for ix in range(sstLons.shape[1]):
                        stateSST[iy,ix]=np.corrcoef(sst_temp[:,0,iy,ix],ST)[0][1]
                fig = plt.figure()
                ax11 = plt.subplot(121);ax12 = plt.subplot(122);
                ax11.text(STyrs[-5],2.25,str(round(STr,3)),fontsize=18,fontweight='bold')
                m1 = Basemap(projection='mill',lon_0=180,lat_0=0,llcrnrlat=-34,urcrnrlat=34,
                 llcrnrlon=60,urcrnrlon=340,ax=ax12); m1.drawcoastlines();m1.drawcountries();
                ax11.plot(STyrs,ensoAll_temp,'k');ax11.plot(STyrs,ST/np.std(ST),'--k',lw=2)
                m1.pcolor(sstLons,sstLats,stateSST,cmap='RdBu_r',shading='flat',latlon=True,vmin=-0.5,vmax=0.5)
                m1.fillcontinents(color='w');m1.drawmeridians(np.arange(60,360,60),linewidth=0,labels=[0,0,0,1])
                m1.drawparallels(np.arange(-30,31,30),linewidth=0,labels=[0,1,0,0])
                fig.set_size_inches(16, 10)
                fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/states/'+dataFile+'_'+STnam+'_wheat_'+str(breakPt)+lag[0]+'lag.png')
                plt.close()
                fig = plt.figure()
                plt.plot(ST_temp,'o')
                fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/trends/'+dataFile+'_'+STnam+'_wheat_yield'+lag[0]+'lag.png')
                plt.close()
                ST = np.append(STyrs.year[:,np.newaxis],ST.T[:,np.newaxis],axis=1)
                np.savetxt(str('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/anomData/'+dataFile+'_'+STnam+'_wheat_yield'+lag[0]+'lag.csv'),ST,delimiter=',')
            if (STr > HighBound)|(STr < LowBound): 
                states1.append(STnam)
                print(str(STnam +' '+ str(sig)+' Corr:'+ str(STr)))
                statesR1.append(STr)
                statesSig1.append(sig)
                REGp1.append(np.float(CNTRYp['2010'][STnam].values[0])/CNTRYp['2010']['Value'].values[0]*100)
        wheat1 = {'state':states1, 'CorrCoeff':statesR1, 'Sig':statesSig1, 'Percent of National Production':REGp1}
        wheatDF1 = pd.DataFrame(data=wheat1)
        pd.DataFrame.to_csv(wheatDF1,
        '/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/yieldCorr/percentAnom/wheat'+country+str(corrThresh)+'.csv')
        print(' ')
        print(' ')

                                            
    #---------------#---------------#---------------#---------------#-------------#
    #---------------------------------- Maize ------------------------------------#
    #---------------#---------------#---------------#---------------#-------------#
    if (crop != 'wheat')&(crop != 'soy'):
        ensoAll = np.array(enso34tab.loc[str(stYr[1]):str(endYr[1]-1)])
        ensoAll = ensoAll[:,0:12]
        if country=='BR':
            maize = BRmaize
        else: maize = 'maize'
        
        #read in the yield data
        CNTRYyld = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/AmENSO_yld/AmENSO_'+country+'_'+maize+'_halfPercMaskGDD2all'+stORreg+'.csv', 
                         header=0)
        CNTRYyld.set_index([CNTRYyld.index.year,'state'],drop=True,inplace=True)    
        CNTRYyld = CNTRYyld['yldAnomGau'];CNTRYyld=CNTRYyld.sort_index()
#        CNTRYyld = CNTRYyld.sort()
        
        #read in the harvested area data
        CNTRYha = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/HarvestedArea/HistoricalData/'+dataFile+'/'+country+'_'+maize+'_'+stORreg+'.csv',
                                      header=0)
        CNTRYha = CNTRYha.sort()
        
        #read in the total production data
        CNTRYp = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/production/HistoricalData/'+dataFile+'/'+country+'_'+maize+'_'+stORreg+'.csv',
                                      header=0)
        CNTRYp = CNTRYp.sort(); #CNTRYp = np.array(CNTRYp[str(yldSt[1]):str(yldEnd[1])]['Value'])
#        if dataFile == 'BR_CONAB':
#            CNTRYp['Value']=CNTRYp['Value']*1000.
        #NaNs coded as 0 in the data
        CNTRYyld[CNTRYyld==0]=np.nan
        CNTRYha[CNTRYha==0]=np.nan
        
        REGp1 = []
        states = []
        statesR1 = []
        statesSig1 = []
        
        #First explore the correlations
        
        #Read in the SST anomaly data
        sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
        #T is months since 1960-01-01, and I want to start at 1950 with one year lag
        sst=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-stYr[1]))&
                                    (sstData.variables['T'][:]<=12*(endYr[1]-1960)),:,
                                    (sstData.variables['Y'][:]<45)&
                                    (sstData.variables['Y'][:]>-45),:]
        #pull only november SST anomalies
        sst=sst[(np.array(range(endYr[1]-stYr[1]))*12)+mon[1],...]
        sstLats=sstData.variables['Y'][(sstData.variables['Y'][:]<45)&
                                    (sstData.variables['Y'][:]>-45)]; sstLons=sstData.variables['X'][:]
        sstLons,sstLats=np.meshgrid(sstLons,sstLats)
        stateSST = np.zeros(np.shape(sstLons))    
        
        #linearly detrend each state yields from 1950-80 and 81-2013
        for STnam in np.unique(CNTRYyld.index.get_level_values(1)):
            surrogateR1 = [] #clear previous surrogate correlations
            if np.array(CNTRYyld.loc[(yldSt[1]):(yldEnd[1]),STnam]).dtype == 'O':
                continue #this skips columns with missing values
            ST_temp = CNTRYyld.loc[(yldSt[1]):(yldEnd[1]),STnam]
            ensoAll_temp = ensoAll[np.isfinite(np.array(ST_temp)),ensoMon[1]]
            ensoTS_temp = ensoTS[np.isfinite(np.array(ST_temp)),1,:] #clip the surrogate TS to dates
            sst_temp = sst[np.isfinite(np.array(ST_temp)),...]
            ST_temp = ST_temp.dropna()
            if len(ST_temp) < endYr[1]-stYr[1]-1:
                continue
            STyrs = ST_temp.index
            ST = ST_temp#np.array(signal.detrend(ST_temp,bp=breakPt)/(ST_temp-signal.detrend(ST_temp,bp=breakPt))) #%yield anom = (yld obs - yld exp)/yld exp
            ST_temp = np.ravel(ST_temp)            
            #ST = np.array(ST_temp[smooth:-smooth]-moving_average(ST_temp))/moving_average(ST_temp)            
            for isu in range(numSurrogates): #correlate against surrogate TS to find significance
                surrogateR1.append(spearmanr(ensoTS_temp[:,isu],ST)[0])
                    #sort and find non parametric alphas
            surrogateR1 = np.sort(surrogateR1)
            LowBound = surrogateR1[numSurrogates*corrThresh//2]
            HighBound = surrogateR1[-numSurrogates*corrThresh//2]                        
            STr = spearmanr(ST,ensoAll_temp)[0]
            sig = spearmanr(ST,ensoAll_temp)[1]
            #ST = np.array(ST_temp)
            #STr = spearmanr(ST[1:]-ST[:-1],ensoAll_temp[1:]-ensoAll_temp[:-1])[0]
            #sig = spearmanr(ST[1:]-ST[:-1],ensoAll_temp[1:]-ensoAll_temp[:-1])[1]
            if printState is 'Y':
                for iy in range(sstLats.shape[0]):
                    for ix in range(sstLons.shape[1]):
                        stateSST[iy,ix]=np.corrcoef(sst_temp[:,0,iy,ix],ST)[0][1]
                fig = plt.figure()
                ax11 = plt.subplot(121);ax12 = plt.subplot(122);
                ax11.text(STyrs[-5],2.25,str(round(STr,3)),fontsize=18,fontweight='bold')
                m1 = Basemap(projection='mill',lon_0=180,lat_0=0,llcrnrlat=-34,urcrnrlat=34,
                 llcrnrlon=60,urcrnrlon=340,ax=ax12); m1.drawcoastlines();m1.drawcountries();
                ax11.plot(STyrs,ensoAll_temp,'k');ax11.plot(STyrs,ST/np.std(ST),'--k',lw=2)
                m1.pcolor(sstLons,sstLats,stateSST,cmap='RdBu_r',shading='flat',latlon=True,vmin=-0.5,vmax=0.5)
                m1.fillcontinents(color='w');m1.drawmeridians(np.arange(60,360,60),linewidth=0,labels=[0,0,0,1])
                m1.drawparallels(np.arange(-30,31,30),linewidth=0,labels=[0,1,0,0])
                fig.set_size_inches(16, 10)
                fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/states/'+dataFile+'_'+STnam+'_'+maize+'_'+str(breakPt)+lag[1]+'lag.png')
                plt.close()
                fig = plt.figure()
                plt.plot(ST_temp,'o')
                fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/trends/'+dataFile+'_'+STnam+'_'+maize+'_yield'+lag[1]+'lag.png')
                plt.close()
                ST = np.append(STyrs.year[:,np.newaxis],ST.T[:,np.newaxis],axis=1)
                np.savetxt(str('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/anomData/'+dataFile+'_'+STnam+'_'+maize+'_yield'+lag[1]+'lag.csv'),ST,delimiter=',')
            if (STr > HighBound)|(STr < LowBound):
                states.append(STnam)        
                print(str(STnam +' '+ str(sig)+' Corr:'+ str(STr)))
                statesR1.append(STr)
                statesSig1.append(sig)
                REGp1.append(np.float(CNTRYp['2010'][STnam].values[0])/CNTRYp['2010']['Value'].values[0]*100)
        maize1 = {'state':states, 'CorrCoeff':statesR1, 'Sig':statesSig1, 'Percent of National Production':REGp1}
        maizeDF1 = pd.DataFrame(data=maize1)
        pd.DataFrame.to_csv(maizeDF1,
        '/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/yieldCorr/percentAnom/'+maize+country+str(corrThresh)+'.csv')

        print(' ')
        print(' ')
     
    
    #---------------#---------------#---------------#---------------#-------------#
    #---------------------------------- Soybean ------------------------------------#
    #---------------#---------------#---------------#---------------#-------------#
    if (crop != 'wheat')&(crop != 'maize'):
        ensoAll = np.array(enso34tab.loc[str(stYr[2]):str(endYr[2]-1)])
        ensoAll = ensoAll[:,0:12]
        
        ensoAll_diff = np.array(enso34tab.loc[str(stYr[2]+1):str(endYr[2])]) #offset by a year for the difference calculation of yield anomalies
        ensoAll_diff = ensoAll_diff[:,0:12]
        
        #read in the yield data
        CNTRYyld = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/AmENSO_yld/AmENSO_'+country+'_soy_halfPercMaskGDD2all'+stORreg+'.csv', 
                         header=0)
        CNTRYyld.set_index([CNTRYyld.index.year,'state'],drop=True,inplace=True)    
        CNTRYyld = CNTRYyld['yldAnomGau'];CNTRYyld=CNTRYyld.sort_index()
#        CNTRYyld = CNTRYyld.sort()
        
        #read in the harvested area data
        CNTRYha = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/HarvestedArea/HistoricalData/'+dataFile+'/'+country+'_soy_'+stORreg+'.csv',
                                      header=0)
        CNTRYha = CNTRYha.sort()
        
        #read in the total production data
        CNTRYp = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/production/HistoricalData/'+dataFile+'/'+country+'_soy_'+stORreg+'.csv',
                                      header=0)
        CNTRYp = CNTRYp.sort(); #CNTRYp = np.array(CNTRYp[str(yldSt[2]):str(yldEnd[2])]['Value'])
#        if dataFile == 'BR_CONAB':
#            CNTRYp=CNTRYp*1000.
        #Read in the SST anomaly data
        sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
        #T is months since 1960-01-01, and I want to start at 1950 with one year lag
        sst=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-stYr[2]))&
                                    (sstData.variables['T'][:]<=12*(endYr[2]-1960)),:,
                                    (sstData.variables['Y'][:]<45)&
                                    (sstData.variables['Y'][:]>-45),:]
        #pull only november SST anomalies
        sst=sst[(np.array(range(endYr[2]-stYr[2]))*12)+mon[2],...]
        
        #repeat for the difference method of yield calculation
        sst_diff=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-stYr[2]-1))&
                                    (sstData.variables['T'][:]<=12*(1+endYr[2]-1960)),:,
                                    (sstData.variables['Y'][:]<45)&
                                    (sstData.variables['Y'][:]>-45),:]
        sst_diff=sst_diff[(np.array(range(endYr[2]-stYr[2]))*12)+mon[2],...]
        
        sstLats=sstData.variables['Y'][(sstData.variables['Y'][:]<45)&
                                    (sstData.variables['Y'][:]>-45)]; sstLons=sstData.variables['X'][:]
        sstLons,sstLats=np.meshgrid(sstLons,sstLats)
        
        #create objects to fill later with different yield anomaly calcs
        stateSST = np.zeros(np.shape(sstLons))
        stateSSTfd = np.zeros(np.shape(sstLons))  
        stateSSTpc = np.zeros(np.shape(sstLons))  
        
        #NaNs coded as 0 in the data
        CNTRYyld[CNTRYyld==0]=np.nan
        CNTRYha[CNTRYha==0]=np.nan
        REGp1 = []
        states3 = []
        statesR1 = []
        statesSig1 = []
        #First explore the correlations
        #linearly detrend each state yields from 1950-80 and 81-2013
        for STnam in np.unique(CNTRYyld.index.get_level_values(1)):
            surrogateR1 = [] #clear previous surrogate correlations
            if (STnam is ' ')|(STnam[:7]=='Unnamed'): continue
            if np.array(CNTRYyld.loc[(yldSt[2]):(yldEnd[2]),STnam]).dtype == 'O':
                continue #this skips columns with missing values
            ST_temp = CNTRYyld.loc[(yldSt[2]):(yldEnd[2]),STnam]
            ensoAll_temp = ensoAll[np.isfinite(np.array(ST_temp)),ensoMon[2]]
            ensoTS_temp = ensoTS[np.isfinite(np.array(ST_temp)),2,:] #clip the surrogate TS to dates
            sst_temp = sst[np.isfinite(np.array(ST_temp)),...]
            ST_temp = ST_temp.dropna()
#            ST_temp_diff= CNTRYyld[str(yldSt[2]):str(yldEnd[2]+1)][STnam]
#            ensoAll_temp_diff = ensoAll_diff[np.isfinite(np.array(ST_temp_diff[:-1])),ensoMon[2]]
#            sst_temp_diff = sst_diff[np.isfinite(np.array(ST_temp_diff[:-1])),...]
#            ST_temp_diff = ST_temp_diff.dropna()
            if len(ST_temp) < endYr[2]-stYr[2]-1:
                continue
            STyrs = ST_temp.index
            ST = ST_temp#np.array(signal.detrend(ST_temp,bp=breakPt)/(ST_temp-signal.detrend(ST_temp,bp=breakPt))) #%yield anom = (yld obs - yld exp)/yld exp
#            STfd = np.array(ST_temp_diff[1:])-np.array(ST_temp_diff[:-1])
#            STpc = (np.array(ST_temp_diff[1:])-np.array(ST_temp_diff[:-1]))/np.array(ST_temp_diff[:-1])
            for isu in range(numSurrogates): #correlate against surrogate TS to find significance
                surrogateR1.append(spearmanr(ensoTS_temp[:,isu],ST)[0])
                    #sort and find non parametric alphas
            surrogateR1 = np.sort(surrogateR1)
            LowBound = surrogateR1[numSurrogates*corrThresh//2]
            HighBound = surrogateR1[-numSurrogates*corrThresh//2]    
#            print(HighBound
            STr = spearmanr(ST,ensoAll_temp)[0]
#            STrfd = spearmanr(STfd,ensoAll_temp_diff)[0]
#            STrpc = spearmanr(STpc,ensoAll_temp_diff)[0]
            
            sig = spearmanr(ST,ensoAll_temp)[1]
#            sigfd = spearmanr(STfd,ensoAll_temp_diff)[1]
#            sigpc = spearmanr(STpc,ensoAll_temp_diff)[1]
            if printState is 'Y':
                for iy in range(sstLats.shape[0]):
                    for ix in range(sstLons.shape[1]):
                        stateSST[iy,ix]=np.corrcoef(sst_temp[:,0,iy,ix],ST)[0][1]
                fig = plt.figure()
                ax11 = plt.subplot(121);ax12 = plt.subplot(122);
                ax11.text(STyrs[-5],2.25,str(round(STr,3)),fontsize=18,fontweight='bold')
                m1 = Basemap(projection='mill',lon_0=180,lat_0=0,llcrnrlat=-34,urcrnrlat=34,
                 llcrnrlon=60,urcrnrlon=340,ax=ax12); m1.drawcoastlines();m1.drawcountries();
                ax11.plot(STyrs,ensoAll_temp,'k');ax11.plot(STyrs,ST/np.std(ST),'--k',lw=2)
                m1.pcolor(sstLons,sstLats,stateSST,cmap='RdBu_r',shading='flat',latlon=True,vmin=-0.5,vmax=0.5)
                m1.fillcontinents(color='w');m1.drawmeridians(np.arange(60,360,60),linewidth=0,labels=[0,0,0,1])
                m1.drawparallels(np.arange(-30,31,30),linewidth=0,labels=[0,1,0,0])
                fig.set_size_inches(16, 10)
                fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/states/'+dataFile+'_'+STnam+'_soy_'+str(breakPt)+lag[2]+'lag.png')
                plt.close()
                fig = plt.figure()
                plt.plot(ST_temp,'o')
                fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/trends/'+dataFile+'_'+STnam+'_soy_yield'+lag[2]+'lag.png')
                plt.close()
                ST = np.append(STyrs.year[:,np.newaxis],ST.T[:,np.newaxis],axis=1)
                np.savetxt(str('/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/anomData/'+dataFile+'_'+STnam+'_soy_yield'+lag[2]+'lag.csv'),ST,delimiter=',')
                
            if (STr > HighBound)|(STr < LowBound):
                states3.append(STnam)
                print(str(STnam +' '+ str(sig)+' Corr:'+ str(STr)))
                statesR1.append(STr)
                statesSig1.append(sig)
                REGp1.append(np.float(CNTRYp['2010'][STnam].values[0])/CNTRYp['2010']['Value'].values[0]*100)
        soy1 = {'state':states3, 'CorrCoeff':statesR1, 'Sig':statesSig1, 'Percent of National Production':REGp1}
        soyDF1 = pd.DataFrame(data=soy1)
        pd.DataFrame.to_csv(soyDF1,
        '/Users/weston/Desktop/Columbia/Research/Results/'+proj+'/Yields/yieldCorr/percentAnom/soy'+country+str(corrThresh)+'.csv')
        print(' ')
        print(' ')
