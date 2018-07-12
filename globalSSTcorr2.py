# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:14:22 2016
@author: weston

Read in the production weighted anomalies, and the percent production anomalies
and correlate them with global SSTs
"""
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
import time
import netCDF4
from scipy import stats
from scipy.stats.stats import pearsonr
from matplotlib import pyplot as plt
plt.ioff() #turn interactive plotting off
start = time.clock()

alpha = .1
numSurrogates = 1000
USwLN = 'LN'
ARsLN = 'LN'
yrMin2 = '2007'; yrMax2 = '2012'#should be a string
#yrMin2 = '1980'; yrMax2 = '1985'
crops = ['maize']#,'soy','maize']#
yMin = -10; yMax = 10
ensoMon = 11 #NOT python indexing
lag = 0
notes = str('November ('+str(lag)+')')
loadNotes = '_HAwgt'
PC = 'PC1'
anomType = 'Production'
dynHA = '' #'' or 'dynHA'

latMax = 25; latMin = -25; lonMax = 275; lonMin = 140

#create surrogate timeseries
#follows from Ebusuzaki (1997) and Schrieber and Shmitz (2000)
#  function call is   ############ surrogate(input TS): ############                  #
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/surrogateData.py').read())
#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#


#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][(sstData.variables['Y'][:]<latMax)&
                            (sstData.variables['Y'][:]>latMin)]; 
sstLons=sstData.variables['X'][(sstData.variables['X'][:]<lonMax)&
                            (sstData.variables['X'][:]>lonMin)]
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

#Read in ENSO information 
enso34tab = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/ENSO/NinoIndex34.csv', 
                     header=1)

insig = dict(linestyle='--',linewidth=1,color='k')
sig = dict(linestyle='-',linewidth=1,color='k')

for crop in crops:
    if (crop is 'wheat')|(crop is 'springwheat')|(crop is 'winterwheat'):
        USw = USwLN
        ARs = ''
    elif (crop is 'soy'):
        USw = ''
        ARs = ARsLN 
    else: 
        USw = ''
        ARs = ''

        
    #Read in the different timeseries
    pcPrDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+anomType+'Anoms/pct'+anomType+dynHA+'_'+crop+'_'+yrMin2+'_'+yrMax2+USw+ARs+loadNotes+'.csv')
    pcPrStDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+anomType+'Anoms/pct'+anomType+dynHA+'ST_'+crop+'_'+yrMin2+'_'+yrMax2+USw+ARs+loadNotes+'.csv')
    pcDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+anomType+'Anoms/PCA_'+crop+'_'+yrMin2+'_'+yrMax2+USw+ARs+'_'+PC+dynHA+loadNotes+'.csv')
    
    sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-1960))&
                            (sstData.variables['T'][:]<=12*((1+2014)-1960)),:,
                            (sstData.variables['Y'][:]<latMax)&(sstData.variables['Y'][:]>latMin),
                            (sstData.variables['X'][:]<lonMax)&(sstData.variables['X'][:]>lonMin)]
    #pull only november SST anomalies
    sstAnom = np.squeeze(sstAnom)
    sstAnom=sstAnom[(np.array(range(2014-1960))*12)+(ensoMon-1),...]

    surrogateR1 = np.zeros([numSurrogates,sstLats.shape[0],sstLons.shape[1]]);
    surrogateR2 = np.zeros([numSurrogates,sstLats.shape[0],sstLons.shape[1]]);
    surrogateR3 = np.zeros([numSurrogates,sstLats.shape[0],sstLons.shape[1]])  #clear previous surrogate correlations
    sstPCcor = np.zeros(np.shape(sstLons)); sstPCp = np.zeros(np.shape(sstLons))
    sstPRcor = np.zeros(np.shape(sstLons)); sstPRp = np.zeros(np.shape(sstLons))
    sstPRSTcor = np.zeros(np.shape(sstLons)); sstPRSTp = np.zeros(np.shape(sstLons))

    for iy in range(sstLats.shape[0]):
        print( iy)
        for ix in range(sstLons.shape[1]):
            tempR1 = []; tempR2 = []; tempR3 = []
            sstPCcor[iy,ix], sstPCp[iy,ix] =pearsonr(sstAnom[(pcDF['year'].min()+lag-1960):(1+pcDF['year'].max()-1960+lag),iy,ix],pcDF['PC1'])
            sstPRcor[iy,ix], sstPRp[iy,ix] =pearsonr(sstAnom[(pcPrDF['year'].min()+lag-1960):(1+pcPrDF['year'].max()-1960+lag),iy,ix],pcPrDF['pctAnom'])
            sstPRSTcor[iy,ix], sstPRSTp[iy,ix] =pearsonr(sstAnom[(pcPrStDF['year'].min()+lag-1960):(1+pcPrStDF['year'].max()-1960+lag),iy,ix],pcPrStDF['pctAnom'])
            #Account for temporal autocorrelation using methods of Ebisuzaki et al.
            for isu in range(numSurrogates): #correlate against surrogate TS to find significance
                ensoTS = surrogate(sstAnom[:,iy,ix])
                ensoTS1 = ensoTS[(pcDF['year'].min()+lag-1960):(1+pcDF['year'].max()-1960+lag)]
                ensoTS2 = ensoTS[(pcPrDF['year'].min()+lag-1960):(1+pcPrDF['year'].max()-1960+lag)]
                ensoTS3 = ensoTS[(pcPrDF['year'].min()+lag-1960):(1+pcPrDF['year'].max()-1960+lag)]
                tempR1.append(pearsonr(ensoTS1,pcDF['PC1'])[0])
                tempR2.append(pearsonr(ensoTS2,pcPrDF['pctAnom'])[0])
                tempR3.append(pearsonr(ensoTS3,pcPrStDF['pctAnom'])[0])
            #sort and find non parametric alphas
            surrogateR1[:,iy,ix] = np.sort(tempR1)
            surrogateR2[:,iy,ix] = np.sort(tempR2)
            surrogateR3[:,iy,ix] = np.sort(tempR3)

    lb1 = surrogateR1[numSurrogates*alpha/2,...]
    hb1 = surrogateR1[-numSurrogates*alpha/2,...]
    lb2 = surrogateR2[numSurrogates*alpha/2,...]
    hb2 = surrogateR2[-numSurrogates*alpha/2,...]
    lb3 = surrogateR3[numSurrogates*alpha/2,...]
    hb3 = surrogateR3[-numSurrogates*alpha/2,...] 

    #Here mask out data that doesn't meet significance levels based on Ebisuzaki methods    
    sstLonsPC = np.copy(sstLons); sstLatsPC = np.copy(sstLats)
    sstLonsPC[(sstPCcor<hb1)&(sstPCcor>lb1)] = np.ma.masked
    sstLatsPC[(sstPCcor<hb1)&(sstPCcor>lb1)] = np.ma.masked  

    sstLonsPR = np.copy(sstLons); sstLatsPR = np.copy(sstLats)
    sstLonsPR[(sstPRcor<hb2)&(sstPRcor>lb2)] = np.ma.masked
    sstLatsPR[(sstPRcor<hb2)&(sstPRcor>lb2)] = np.ma.masked 

    sstLonsPRST = np.copy(sstLons); sstLatsPRST = np.copy(sstLats)
    sstLonsPRST[(sstPRSTcor<hb3)&(sstPRSTcor>lb3)] = np.ma.masked
    sstLatsPRST[(sstPRSTcor<hb3)&(sstPRSTcor>lb3)] = np.ma.masked 

    #Now use the methods of Wilks (2016) to account for multiple hypothesis testing
    #first find the number and location of point you're testing
    #As noted in Wilks (2016), you can account for strong spatial correlations by
    # approximating alpha_FDR = alpha_global*2
    N = sstLats.size
    
    p_FDR_PC = np.copy(np.ravel(sstPCp))
    p_FDR_PC.sort();p_FDR_PC[p_FDR_PC==1.0]=np.nan
    
    p_FDR_PR = np.copy(np.ravel(sstPRp))
    p_FDR_PR.sort();p_FDR_PR[p_FDR_PR==1.0]=np.nan

    p_FDR_PRST = np.copy(np.ravel(sstPRSTp))
    p_FDR_PRST.sort();p_FDR_PRST[p_FDR_PRST==1.0]=np.nan
    
    for i in range(N):
        if (p_FDR_PC[i]>((i+1)/np.float(N)*alpha*2)): 
            p_FDR_PC[i]=0        
        if (p_FDR_PR[i]>((i+1)/np.float(N))*alpha*2): 
            p_FDR_PR[i]=0
        if (p_FDR_PRST[i]>((i+1)/np.float(N))*alpha*2): 
            p_FDR_PRST[i]=0

    #Determine the p_FDR values        
    p_FDR_PC = np.nanmax(p_FDR_PC)
    p_FDR_PR = np.nanmax(p_FDR_PR)
    p_FDR_PRST = np.nanmax(p_FDR_PRST)

    #apply the values
    sstLonsPC[sstPCp>p_FDR_PC] = np.ma.masked
    sstLatsPC[sstPCp>p_FDR_PC] = np.ma.masked 
    sstLonsPR[sstPRp>p_FDR_PR] = np.ma.masked
    sstLatsPR[sstPRp>p_FDR_PR] = np.ma.masked 
    sstLonsPRST[sstPRSTp>p_FDR_PRST] = np.ma.masked
    sstLatsPRST[sstPRSTp>p_FDR_PRST] = np.ma.masked 

    pcEN_1 = np.array(pcDF.loc[pcDF['EN']==-1]['PC1'])    
    pcEN0 = np.array(pcDF.loc[pcDF['EN']==0]['PC1']) 
    pcEN1 = np.array(pcDF.loc[pcDF['EN']==1]['PC1'])
    pcEN_1sig = stats.ttest_1samp(pcEN_1,0).pvalue
    pcEN0sig = stats.ttest_1samp(pcEN0,0).pvalue
    pcEN1sig = stats.ttest_1samp(pcEN1,0).pvalue
    pcLN_1 = np.array(pcDF.loc[pcDF['LN']==-1]['PC1'])    
    pcLN0 = np.array(pcDF.loc[pcDF['LN']==0]['PC1']) 
    pcLN1 = np.array(pcDF.loc[pcDF['LN']==1]['PC1'])
    pcLN_1sig = stats.ttest_1samp(pcLN_1,0).pvalue
    pcLN0sig = stats.ttest_1samp(pcLN0,0).pvalue
    pcLN1sig = stats.ttest_1samp(pcLN1,0).pvalue
    pcEN = [pcEN_1,pcEN0,pcEN1];pcLN = [pcLN_1,pcLN0,pcLN1]
        
    pcPrEN_1 = np.array(pcPrDF.loc[pcPrDF['EN']==-1]['pctAnom'])    
    pcPrEN0 = np.array(pcPrDF.loc[pcPrDF['EN']==0]['pctAnom']) 
    pcPrEN1 = np.array(pcPrDF.loc[pcPrDF['EN']==1]['pctAnom'])
    pcPrEN_1sig = stats.ttest_1samp(pcPrEN_1,0).pvalue
    pcPrEN0sig = stats.ttest_1samp(pcPrEN0,0).pvalue
    pcPrEN1sig = stats.ttest_1samp(pcPrEN1,0).pvalue
    pcPrLN_1 = np.array(pcPrDF.loc[pcPrDF['LN']==-1]['pctAnom'])    
    pcPrLN0 = np.array(pcPrDF.loc[pcPrDF['LN']==0]['pctAnom']) 
    pcPrLN1 = np.array(pcPrDF.loc[pcPrDF['LN']==1]['pctAnom'])
    pcPrLN_1sig = stats.ttest_1samp(pcPrLN_1,0).pvalue
    pcPrLN0sig = stats.ttest_1samp(pcPrLN0,0).pvalue
    pcPrLN1sig = stats.ttest_1samp(pcPrLN1,0).pvalue
    pcPrEN = [pcPrEN_1,pcPrEN0,pcPrEN1];pcPrLN = [pcPrLN_1,pcPrLN0,pcPrLN1]    
    
    pcPrStEN_1 = np.array(pcPrStDF.loc[pcPrStDF['EN']==-1]['pctAnom'])  
    pcPrStEN0 = np.array(pcPrStDF.loc[pcPrStDF['EN']==0]['pctAnom'])
    pcPrStEN1 = np.array(pcPrStDF.loc[pcPrStDF['EN']==1]['pctAnom'])
    pcPrStEN_1sig = stats.ttest_1samp(pcPrStEN_1,0).pvalue
    pcPrStEN0sig = stats.ttest_1samp(pcPrStEN0,0).pvalue
    pcPrStEN1sig = stats.ttest_1samp(pcPrStEN1,0).pvalue
    pcPrStLN_1 = np.array(pcPrStDF.loc[pcPrStDF['LN']==-1]['pctAnom']) 
    pcPrStLN0 = np.array(pcPrStDF.loc[pcPrStDF['LN']==0]['pctAnom'])
    pcPrStLN1 = np.array(pcPrStDF.loc[pcPrStDF['LN']==1]['pctAnom'])
    pcPrStLN_1sig = stats.ttest_1samp(pcPrStLN_1,0).pvalue
    pcPrStLN0sig = stats.ttest_1samp(pcPrStLN0,0).pvalue
    pcPrStLN1sig = stats.ttest_1samp(pcPrStLN1,0).pvalue
    pcPrStEN = [pcPrStEN_1,pcPrStEN0,pcPrStEN1];pcPrStLN = [pcPrStLN_1,pcPrStLN0,pcPrStLN1]
    pcPrStEN_sig = [pcPrStEN_1sig,pcPrStEN0sig,pcPrStEN1sig] ;pcPrStLN_sig = [pcPrStLN_1sig,pcPrStLN0sig,pcPrStLN1sig]    
    pcPrEN_sig = [pcPrEN_1sig,pcPrEN0sig,pcPrEN1sig] ;pcPrLN_sig = [pcPrLN_1sig,pcPrLN0sig,pcPrLN1sig]    
    pcEN_sig = [pcEN_1sig,pcEN0sig,pcEN1sig] ;pcLN_sig = [pcLN_1sig,pcLN0sig,pcLN1sig]    
    
    #convert significances into linestyles
    for il in range(3):
        if pcEN_sig[il]>alpha: pcEN_sig[il]=':'
        else: pcEN_sig[il]='-'
        if pcLN_sig[il]>alpha: pcLN_sig[il]=':'
        else: pcLN_sig[il]='-'

        if pcPrEN_sig[il]>alpha: pcPrEN_sig[il]=':'
        else: pcPrEN_sig[il]='-'
        if pcPrLN_sig[il]>alpha: pcPrLN_sig[il]=':'
        else: pcPrLN_sig[il]='-'

        if pcPrStEN_sig[il]>alpha: pcPrStEN_sig[il]=':'
        else: pcPrStEN_sig[il]='-'
        if pcPrStLN_sig[il]>alpha: pcPrStLN_sig[il]=':'
        else: pcPrStLN_sig[il]='-'
    
    fig = plt.figure()
    ax11 = plt.subplot(331);ax12 = plt.subplot(332);ax13 = plt.subplot(333);
    ax21 = plt.subplot(334);ax22 = plt.subplot(335);ax23 = plt.subplot(336);
    ax31 = plt.subplot(337);ax32 = plt.subplot(338);ax33 = plt.subplot(339);
    ax11.set_title(u"El Niño life-cycle\n EN -1        EN 0        EN 1",size=18)
    ax12.set_title(u"La Niña life-cycle\n LN -1        LN 0        LN 1",size=18)
    ax13.set_title(u"Correlation with\n SST anomalies",size=18)
    ax13.yaxis.labelpad = 15;ax23.yaxis.labelpad = 15;ax33.yaxis.labelpad = 15;ax21.yaxis.labelpad = 25
    y1 = ax13.set_ylabel(str('All states \n'+notes),size=16);y1.set_rotation(90);ax13.yaxis.set_label_position('left')    
    y2 = ax23.set_ylabel(str('Correlated states \n'+notes),size=16);y2.set_rotation(90);ax23.yaxis.set_label_position('left')  
    y3 = ax33.set_ylabel(str('PC '+PC[2]+'\n'+notes),size=16);y3.set_rotation(90);ax33.yaxis.set_label_position('left')  
    y4 = ax21.set_ylabel('Percent production anomaly',size=16);
    
    
    m1 = Basemap(projection='merc',llcrnrlat=latMin,llcrnrlon=lonMin,urcrnrlat=latMax,urcrnrlon=lonMax,ax=ax13); m1.drawcoastlines();m1.drawcountries();
    m2 = Basemap(projection='merc',llcrnrlat=latMin,llcrnrlon=lonMin,urcrnrlat=latMax,urcrnrlon=lonMax,ax=ax23); m2.drawcoastlines();m2.drawcountries();
    m3 = Basemap(projection='merc',llcrnrlat=latMin,llcrnrlon=lonMin,urcrnrlat=latMax,urcrnrlon=lonMax,ax=ax33); m3.drawcoastlines();m3.drawcountries();
    
    
    im1=m1.pcolor(sstLons,sstLats,sstPRcor,cmap='RdBu_r',shading='flat',latlon=True,vmin=-.8,vmax=.8)
    im2=m2.pcolor(sstLons,sstLats,sstPRSTcor,cmap='RdBu_r',shading='flat',latlon=True,vmin=-.8,vmax=.8)
    im3 = m3.pcolor(sstLons,sstLats,sstPCcor,cmap='RdBu_r',shading='flat',latlon=True,vmin=-.8,vmax=.8)
    lonPR, latPR = m1(sstLonsPR,sstLatsPR); m1.plot(lonPR,latPR,'ko',markersize=.5)
    lonPRST, latPRST = m2(sstLonsPRST,sstLatsPRST); m2.plot(lonPRST,latPRST,'ko',markersize=.5)
    lonPC, latPC = m3(sstLonsPC,sstLatsPC); m3.plot(lonPC,latPC,'ko',markersize=.5)

    ax11.set_ylim([yMin,yMax]);ax12.set_ylim([yMin,yMax])
    ax21.set_ylim([yMin,yMax]);ax22.set_ylim([yMin,yMax])
    ax31.set_ylim([yMin,yMax]);ax32.set_ylim([yMin,yMax])
    
    box11 = ax11.boxplot(pcPrEN,showmeans=True, patch_artist=True);
    box12 = ax12.boxplot(pcPrLN,showmeans=True, patch_artist=True)
    box21 = ax21.boxplot(pcPrStEN,showmeans=True, patch_artist=True)
    box22 = ax22.boxplot(pcPrStLN,showmeans=True, patch_artist=True)
    box31 = ax31.boxplot(pcEN,showmeans=True,patch_artist=True)   
    box32 = ax32.boxplot(pcLN,showmeans=True,patch_artist=True)
    
    for ib in range(3):
        box11['boxes'][ib].set(color='#000000',ls=pcPrEN_sig[ib])
        box11['boxes'][ib].set(facecolor='#FFFFFF')
        box11['whiskers'][ib].set(color='#000000')
        box12['boxes'][ib].set(color='#000000',ls=pcPrLN_sig[ib])
        box12['boxes'][ib].set(facecolor='#FFFFFF')
        box12['whiskers'][ib].set(color='#000000')    
        
        box21['boxes'][ib].set(color='#000000',ls=pcPrStEN_sig[ib])
        box21['boxes'][ib].set(facecolor='#FFFFFF')
        box21['whiskers'][ib].set(color='#000000')
        box22['boxes'][ib].set(color='#000000',ls=pcPrStLN_sig[ib])
        box22['boxes'][ib].set(facecolor='#FFFFFF')
        box22['whiskers'][ib].set(color='#000000')  
        
        box31['boxes'][ib].set(color='#000000',ls=pcEN_sig[ib])
        box31['boxes'][ib].set(facecolor='#FFFFFF')
        box31['whiskers'][ib].set(color='#000000')
        box32['boxes'][ib].set(color='#000000',ls=pcLN_sig[ib])
        box32['boxes'][ib].set(facecolor='#FFFFFF')
        box32['whiskers'][ib].set(color='#000000')  
        
    ax11.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k')
    ax12.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k')
    ax21.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k')
    ax22.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k')
    ax31.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k')
    ax32.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k')
    m1.fillcontinents(color='w');m1.drawmeridians(np.arange(60,360,60),linewidth=0,labels=[0,0,0,1])
    m2.fillcontinents(color='w');m2.drawmeridians(np.arange(60,360,60),linewidth=0,labels=[0,0,0,1])
    m3.fillcontinents(color='w');m3.drawmeridians(np.arange(60,360,60),linewidth=0,labels=[0,0,0,1])
    fig.set_size_inches(16, 10)
    ax11.xaxis.set_ticks([]);ax21.xaxis.set_ticks([]);ax31.xaxis.set_ticks([])    
    ax12.xaxis.set_ticks([]);ax22.xaxis.set_ticks([]);ax32.xaxis.set_ticks([])
    cb33 = m3.colorbar(im3,"bottom", size="10%", pad="20%");
    
    fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/AmENSO_yields/Yields/yldLifeCycle/'+anomType+'/'+crop+'_'+yrMin2+yrMax2+'lifeCycle'+ARs+USw+notes+PC+dynHA+loadNotes+'.png')
    plt.close()
