# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:14:22 2016
@author: westonanderson

Read in the production weighted anomalies, and the percent production anomalies
and correlate them with global SSTs
"""
import numpy as np
import pandas as pd
import time
from scipy import stats
from matplotlib import pyplot as plt
start = time.clock()
plt.ioff() #turn interactive plotting off

corrThresh = .1
USwLN = 'LN'
ARsLN = 'LN'
yrMin = '2007' #should be a string
yrMax = '2012'#should be a string
crops = ['springwheat']#,'springwheat']
yMin = -15; yMax = 15
notes = '_HAwgt'
PC = 'PC1'
PCsign = 1 #1 or -1
cnt = 'cnt'

#create surrogate timeseries
#follows from Ebusuzaki (1997) and Schrieber and Shmitz (2000)
def surrogate(ai):
    ak = np.fft.fft(ai) #DFFT
    k = np.size(ai)
    rk = [None]*k
    rk[0] = ak[0]*np.sin(ak[0]) #amp zero frequency (with sign)
    if k%2==0: #if even
        for i in range(1,k/2):
            thetak = np.pi*2.*np.random.uniform(0,1)   #generate random phase      
            rk[i] = np.abs(ak[i])*np.exp(thetak*1j)    #apply random phase to positive freq
            rk[-(i)] = np.abs(ak[-(i+1)])*np.exp(thetak*-1j) #apply random phase to negative freq
        thetak = np.pi*2.*np.random.uniform(0,1) #generate random phase    
        rk[k/2]=np.sqrt(2.)*np.abs(ak[k/2])*np.cos(thetak) #amp the nyquist freq
    else: #if odd
        for i in range(1,(k+1)/2):
            thetak = np.pi*2.*np.random.uniform(0,1)      #generate random phase       
            rk[i] = np.abs(ak[i])*np.exp(thetak*1j) #apply random phase to positive freq
            rk[-(i)] = np.abs(ak[i])*np.exp(thetak*-1j) #apply random phase to negative freq
    rj = np.fft.ifft(rk) #inverse fft
    return rj.real

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
    YldPrDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/ProductionAnoms/PCA_'+crop+'YldAll'+cnt+'_'+yrMin+'_'+yrMax+USw+ARs+notes+'.csv')
    YldPrStDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/ProductionAnoms/PCA_'+crop+'YldST'+cnt+'_'+yrMin+'_'+yrMax+USw+ARs+notes+'.csv')
    YldDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/ProductionAnoms/PCA_'+crop+'Yld'+cnt+'_'+yrMin+'_'+yrMax+USw+ARs+'_'+PC+notes+'DynHA.csv')         
    YldPrDF['pctAnom']=YldPrDF['pctAnom']*100    
    YldPrStDF['pctAnom']=YldPrStDF['pctAnom']*100  
    YldDF['pctAnom']=YldDF['pctAnom']*PCsign*100
    
    fig = plt.figure()
    ax11 = plt.subplot(321);ax12 = plt.subplot(322);
    ax21 = plt.subplot(323);ax22 = plt.subplot(324);
    ax31 = plt.subplot(325);ax32 = plt.subplot(326);
      
    harvs = np.unique(YldPrDF['harvest'])
    cntry = np.unique(YldPrDF['country'])
    i=0
    xTicks = []
    xLabs = []
    lstPos=[]
    for iEN in range(-1,2,1):
        for iHar in harvs:
            i=i+1
            for icntry in cntry:
                iFlw = np.unique(YldDF[YldDF['country']==icntry]['flwrSeas'])[0]
                imo=iEN*12 + iFlw
                YldPrStEN = np.array(YldPrDF.loc[YldPrDF['EN']==iEN].loc[YldPrDF['harvest']==iHar].loc[YldPrDF['country']==icntry]['pctAnom']) 
                YldPrStLN = np.array(YldPrDF.loc[YldPrDF['LN']==iEN].loc[YldPrDF['harvest']==iHar].loc[YldPrDF['country']==icntry]['pctAnom']) 
                YldPrStLN=YldPrStLN[~np.isnan(YldPrStLN)]; YldPrStEN=YldPrStEN[~np.isnan(YldPrStEN)]            
                if np.size(YldPrStEN)==0: continue
                pos = imo
                if pos == lstPos:
                    pos=pos+1
                    lstPos=pos
                else: lstPos=pos
                YldPrStENsig = stats.ttest_1samp(YldPrStEN,0).pvalue
                YldPrStLNsig = stats.ttest_1samp(YldPrStLN,0).pvalue
                if YldPrStENsig >corrThresh: YldPrStENsig=':'
                else: YldPrStENsig ='-'
                box11 = ax11.boxplot(YldPrStEN,showmeans=True, positions=[pos],widths=0.8); #[(iEN*12)+iHar]
                box11['boxes'][0].set(color='#000000',ls=YldPrStENsig)
                box11['whiskers'][0].set(color='#000000')
                box11['whiskers'][1].set(color='#000000')
                if YldPrStLNsig >corrThresh: YldPrStLNsig=':'
                else: YldPrStLNsig ='-'
                box12 = ax12.boxplot(YldPrStLN,showmeans=True, patch_artist=True,positions=[pos],widths=.8);
                box12['boxes'][0].set(color='#000000',ls=YldPrStLNsig)
                box12['boxes'][0].set(facecolor='#FFFFFF')
                box12['whiskers'][0].set(color='#000000')
                box12['whiskers'][1].set(color='#000000')
                if icntry=='United States': cntAbr = 'US'
                elif icntry=='Argentina': cntAbr = 'AR'
                elif icntry=='Brazil': cntAbr = 'BR'
                elif icntry=='Canada': cntAbr = 'CA'
                elif icntry=='Mexico': cntAbr = 'MX'
                lab = str(cntAbr+' ('+str(iHar)+')'+'('+str(iFlw)+')')
                xTicks = np.append(xTicks,pos)
                xLabs = np.append(xLabs,lab)
    ax11.set_xlim([xTicks.min()-1,xTicks.max()+1]);ax12.set_xlim([xTicks.min()-1,xTicks.max()+1])
    ax11.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k');ax12.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k');
#    ax11.axvline(x=1.5);ax11.axvline(x=13.5);ax11.axvline(x=25.5);
#    ax12.axvline(x=1.5);ax12.axvline(x=13.5);ax12.axvline(x=25.5) 
    ax11.set_ylim([-25,25]);ax12.set_ylim([-25,25])


    harvs = np.unique(YldPrStDF['harvest'])
    i=0
    lstPos=[]
    for iEN in range(-1,2,1):
        for iHar in harvs:
            i=i+1
            for icntry in cntry:
                iFlw = np.unique(YldDF[YldDF['country']==icntry]['flwrSeas'])[0]
                imo=iEN*12 + iFlw                
                YldPrStEN = np.array(YldPrStDF.loc[YldPrStDF['EN']==iEN].loc[YldPrStDF['harvest']==iHar].loc[YldPrStDF['country']==icntry]['pctAnom']) 
                YldPrStLN = np.array(YldPrStDF.loc[YldPrStDF['LN']==iEN].loc[YldPrStDF['harvest']==iHar].loc[YldPrStDF['country']==icntry]['pctAnom']) 
                YldPrStLN=YldPrStLN[~np.isnan(YldPrStLN)]; YldPrStEN=YldPrStEN[~np.isnan(YldPrStEN)]            
                if np.size(YldPrStEN)==0: continue
                pos = imo
                if pos == lstPos:
                    pos=pos+1
                    lstPos=pos
                else: lstPos=pos          
                YldPrStENsig = stats.ttest_1samp(YldPrStEN,0).pvalue
                YldPrStLNsig = stats.ttest_1samp(YldPrStLN,0).pvalue
     #           ax31.violinplot(YldPrStEN, positions=[i])
     #           ax32.violinplot(YldPrStLN, positions=[i])
                if YldPrStENsig >corrThresh: YldPrStENsig=':'
                else: YldPrStENsig ='-'
                box21 = ax21.boxplot(YldPrStEN,showmeans=True, positions=[pos],widths=0.8); #[(iEN*12)+iHar]
                box21['boxes'][0].set(color='#000000',ls=YldPrStENsig)
    #            box21['boxes'][0].set(facecolor='#FFFFFF')
                box21['whiskers'][0].set(color='#000000')
                box21['whiskers'][1].set(color='#000000') 
                if YldPrStLNsig >corrThresh: YldPrStLNsig=':'
                else: YldPrStLNsig ='-'
                box22 = ax22.boxplot(YldPrStLN,showmeans=True, patch_artist=True,positions=[pos],widths=.8);
                box22['boxes'][0].set(color='#000000',ls=YldPrStLNsig)
                box22['boxes'][0].set(facecolor='#FFFFFF')
                box22['whiskers'][0].set(color='#000000')
                box22['whiskers'][1].set(color='#000000')
    ax21.set_xlim([xTicks.min()-1,xTicks.max()+1]);ax22.set_xlim([xTicks.min()-1,xTicks.max()+1])
    ax21.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k');ax22.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k');
#    ax21.axvline(x=1.5);ax21.axvline(x=13.5);ax21.axvline(x=25.5);
#    ax22.axvline(x=1.5);ax22.axvline(x=13.5);ax22.axvline(x=25.5) 
    ax21.set_ylim([-25,25]);ax22.set_ylim([-25,25])

    harvs = np.unique(YldDF['harvest'])
    i=0
    lstPos=[]
    for iEN in range(-1,2,1):
        for iHar in harvs:
            i=i+1
            for icntry in cntry:
                iFlw = np.unique(YldDF[YldDF['country']==icntry]['flwrSeas'])[0]
                imo=iEN*12 + iFlw                
                YldPrStEN = np.array(YldDF.loc[YldDF['EN']==iEN].loc[YldDF['harvest']==iHar].loc[YldDF['country']==icntry]['pctAnom']) 
                YldPrStLN = np.array(YldDF.loc[YldDF['LN']==iEN].loc[YldDF['harvest']==iHar].loc[YldDF['country']==icntry]['pctAnom']) 
                YldPrStLN=YldPrStLN[~np.isnan(YldPrStLN)]; YldPrStEN=YldPrStEN[~np.isnan(YldPrStEN)]            
                if np.size(YldPrStEN)==0: continue
                pos = imo
                if pos == lstPos:
                    pos=pos+1
                    lstPos=pos
                else: lstPos=pos            
                YldPrStENsig = stats.ttest_1samp(YldPrStEN,0).pvalue
                YldPrStLNsig = stats.ttest_1samp(YldPrStLN,0).pvalue
     #           ax31.violinplot(YldPrStEN, positions=[i])
     #           ax32.violinplot(YldPrStLN, positions=[i])
                if YldPrStENsig >corrThresh: YldPrStENsig=':'
                else: YldPrStENsig ='-'
                box31 = ax31.boxplot(YldPrStEN,showmeans=True, positions=[pos],widths=0.8); #[(iEN*12)+iHar]
                box31['boxes'][0].set(color='#000000',ls=YldPrStENsig)
    #            box31['boxes'][0].set(facecolor='#FFFFFF')
                box31['whiskers'][0].set(color='#000000')  
                box31['whiskers'][1].set(color='#000000')
                if YldPrStLNsig >corrThresh: YldPrStLNsig=':'
                else: YldPrStLNsig ='-'
                box32 = ax32.boxplot(YldPrStLN,showmeans=True, patch_artist=True,positions=[pos],widths=.8);
                box32['boxes'][0].set(color='#000000',ls=YldPrStLNsig)
                box32['boxes'][0].set(facecolor='#FFFFFF')
                box32['whiskers'][0].set(color='#000000')
                box32['whiskers'][1].set(color='#000000')

    ax31.set_xlim([xTicks.min()-1,xTicks.max()+1]);ax32.set_xlim([xTicks.min()-1,xTicks.max()+1])
    ax31.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k');ax32.axhline(y=0,xmin=0,xmax=4,ls='-.',color='k');
#    ax31.axvline(x=1.5);ax31.axvline(x=13.5);ax31.axvline(x=25.5);
#    ax32.axvline(x=1.5);ax32.axvline(x=13.5);ax32.axvline(x=25.5)  
    ax31.set_ylim([-25,25]);ax32.set_ylim([-25,25])
    ax31.xaxis.set_ticks(xTicks)
    ax31.xaxis.set_ticklabels(xLabs,rotation=90,fontsize=12)
    ax32.xaxis.set_ticks(xTicks)
    ax32.xaxis.set_ticklabels(xLabs,rotation=90,fontsize=12)
    ax11.xaxis.set_ticks([]);ax12.xaxis.set_ticks([])
    ax21.xaxis.set_ticks([]);ax22.xaxis.set_ticks([])  

    ax11.set_title(u"El Niño life-cycle\n EN -1                 EN 0                 EN 1",size=18)
    ax12.set_title(u"La Niña life-cycle\n LN -1                 LN 0                 LN 1",size=18)

    ax12.yaxis.labelpad = 25;ax22.yaxis.labelpad = 25;ax32.yaxis.labelpad = 25;ax21.yaxis.labelpad = 25
    y1 = ax12.set_ylabel(str('All states'),size=16);y1.set_rotation(-90);ax12.yaxis.set_label_position('right')    
    y2 = ax22.set_ylabel(str('Correlated states'),size=16);y2.set_rotation(-90);ax22.yaxis.set_label_position('right')  
    y3 = ax32.set_ylabel(str('EOF '+PC[2]),size=16);y3.set_rotation(-90);ax32.yaxis.set_label_position('right')  
    y4 = ax21.set_ylabel('Percent yield anomaly',size=18);y3.set_rotation(-90);  
    
    
    fig.set_size_inches(16, 10)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/Yields/yldLifeCycle/Production/'+crop+'_'+yrMin+'_'+yrMax+'_yldLifeCycle'+cnt+ARs+USw+notes+PC+'.pdf')
    plt.close()
