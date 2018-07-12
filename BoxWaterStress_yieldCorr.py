# -*- coding: utf-8 -*-
"""
Created on Mon May 23 10:42:05 2016

@author: westonanderson
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

months = [-2,-1] #previous Nov
state = 'MT'
SWC = 'cerrado'
lag = ''
anom = ''
cropNam = 'soy'
yldAnom = 'yldAnom'
years = [1950,2011]
rustYrs = []#'2002','2003','2006'] #2005 also?

swcAll = np.genfromtxt('/Volumes/Data_Archive/Results/soilWaterBalance/'+SWC+'TS'+state+anom+lag+str(years[0])+'_'+str(years[1])+'.txt')

crop = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/AmENSO_yld/AmENSO_BR_'+cropNam+'_allGDD2allstate.csv')
crop = crop.set_index('state',drop=True,append=True)
yields = crop[yldAnom].loc[:,state]

swc = np.zeros([np.shape(months)[0],years[1]-years[0]-1])*np.nan

for im in range(np.shape(months)[0]):
    mon = months[im]
    if mon<0:
        yields = yields[str(years[0]+1):str(years[1]-1)]
        index = np.array(range((swcAll.shape[0])/12))*12 + mon
        index = index[1:]
    else:
        yields = yields[str(years[0]):str(years[1]-1)]
        index = np.array(range((swcAll.shape[0])/12))*12 + mon
    swc[im,...] = swcAll[index]
swc = np.nanmean(swc,axis=0)

for yr in rustYrs:
    yields[yr] = np.nan
nonNan = np.array(~np.isnan(yields))
swc = swc[nonNan]
yields = yields[nonNan]

print np.corrcoef(swc,yields)[0][1]


