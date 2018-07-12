# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:09:34 2015

@author: westonanderson
"""
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap

start = time.clock()

month = 10 #note that this is in python indexing
stYr=1900 #1860 #(1860 full data)
endYr=2014

#Read in ENSO information for later
enso34tab = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/ENSO/Nino34_Extended.csv', 
                     header=1)
ensoAll = np.array(enso34tab.loc[str(stYr):str(endYr-1)])
ensoAll =  np.reshape(ensoAll,[np.size(ensoAll)/12,12])
ensoThresh = -np.std(ensoAll[:,month])/2.
#ensoThresh = -0.48
coldEvents = np.ones([1,36])*np.nan
coldYears = np.zeros(1,dtype=int)
coldEvents1 = np.ones([1,36])*np.nan
coldYears1 = np.zeros(1,dtype=int)
lstYr = 0
for row in range(1,ensoAll.shape[0]):
        if ensoAll[row,month]<=ensoThresh: #mark the threshold variables
            if lstYr == row-1: #if last year was a cold year skip this one, since it will be included as the t+1 
                if ensoAll[lstYr,month]<=ensoThresh:
                    lstYr=row
                    coldYears1 = np.append(coldYears,[stYr+row],axis=0)                    
                    continue
                else: continue
            event = np.ravel(ensoAll[row-1:row+2,:]) 
            coldEvents = np.append(coldEvents,[event],axis=0) 
            coldYears = np.append(coldYears,[stYr+row],axis=0) #cold years are the ENSO years, not the t-1 yrs 
            lstYr = row
coldEvents = coldEvents[1:,:] #strip the empty first event
coldYears = coldYears[1:] #strip the empty first event    
    

#Read in ENSO information for later
enso34tab = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/ENSO/Nino34_Extended.csv', 
                     header=1)
ensoAll = np.array(enso34tab.loc[str(stYr):str(endYr-1)])
ensoAll =  np.reshape(ensoAll,[np.size(ensoAll)/12,12])
ensoThresh = np.std(ensoAll[:,month])/2.
#ensoThresh = 0.48
warmEvents = np.ones([1,36])*np.nan
warmYears = np.zeros(1,dtype=int)
lstYr = 0
for row in range(1,ensoAll.shape[0]):
        if ensoAll[row,month]>=ensoThresh: #mark the threshold variables 
            if lstYr == row-1: #if last year was a warm year skip this one, since it will be included as the t+1 
                if ensoAll[lstYr,month]>=ensoThresh:
                    lstYr=row;continue
                else: continue 
            event = np.ravel(ensoAll[row-1:row+2,:]) 
            warmEvents = np.append(warmEvents,[event],axis=0) 
            warmYears = np.append(warmYears,[stYr+row],axis=0) #warm years are the ENSO years, not the t-1 yrs 
            lstYr = row 
warmEvents = warmEvents[1:,:] #strip the empty first event
warmYears = warmYears[1:] #strip the empty first event


#############**##################**####################**#############**##################**####################
#############**##################**####################**#############**##################**####################

fig = plt.figure(); 
ax1=plt.subplot(212);ax2=plt.subplot(211)

#plot the ensemble of enso evolution in the top panel

ax1.plot(coldEvents.T, ls='--', color='grey')
ax1.set_xticks(np.arange(5,36,6))
ax1.set_xticklabels(['Jul -1','Jan 0', 'Jul 0','Jan +1', 'Jul +1',''])
ax1.hlines(0,0,35);ax1.set_xlim(0,35)
ax1.vlines([11,23],-3,3,linewidth=2)
for iy in range(coldYears.size):
    ax1.text(iy*2+2, 3.1,str(coldYears[iy]),
         horizontalalignment='right',color='grey')
ax1.plot(np.mean(coldEvents,axis=0).T,linewidth=3,color='darkblue')
ax1.text(4,2.6,'LN -1',fontsize=12)#, fontweight='bold')
ax1.text(16,2.6,'LN 0',fontsize=12)#, fontweight='bold')
ax1.text(28,2.6,'LN +1',fontsize=12)#, fontweight='bold')
ax1.set_ylabel('Oceanic Nino Index', size=12)

ax2.plot(warmEvents.T, ls='--', color='grey')
ax2.set_xticks(np.arange(5,36,6))
ax2.set_xticklabels(['Jul -1','Jan 0', 'Jul 0','Jan +1', 'Jul +1',''])
ax2.hlines(0,0,35);ax1.set_xlim(0,35)
ax2.vlines([11,23],-3,3,linewidth=2)
for iy in range(warmYears.size):
    ax2.text(iy*2+2, 3.1,str(warmYears[iy]),
         horizontalalignment='right',color='grey')
         
ax2.plot(np.mean(warmEvents,axis=0).T,linewidth=3,color='darkred')
ax2.text(4,2.6,'EN -1',fontsize=12)#, fontweight='bold')
ax2.text(16,2.6,'EN 0',fontsize=12)#, fontweight='bold')
ax2.text(28,2.6,'EN +1',fontsize=12)#, fontweight='bold')
ax2.set_ylabel('Oceanic Nino Index', size=12)

fig.set_size_inches(10, 10)
plt.show()
fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO/figures/eventSet_long.png')
plt.close()    
#############**##################**####################**#############**##################**####################
#############**##################**####################**#############**##################**####################
print ''
print 'probability of EN in any year: '+ np.str(np.float(np.shape(warmEvents)[0])/(endYr-stYr))
print 'probability of LN in any year: '+ np.str((np.float(sum(coldEvents[:,33]<0.48))+np.float(np.shape(coldEvents)[0]))/(endYr-stYr))
print 'probability of EN -> LN:  ' + np.str(np.float(sum(warmEvents[:,33]<0.48))/np.shape(warmEvents)[0])
print 'probability of LN0 -> LN1:  ' + np.str(np.float(sum(coldEvents[:,33]<0.48))/np.shape(coldEvents)[0])