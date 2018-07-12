# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 13:20:47 2015
@author: westonanderson

Construct regression models using state-wise data and aggregated anomalies.
The mean of each anomaly (including sum precip) is used across the state to match yields

precipitation from the winter season is used directly rather than linearly
relating it first to soil moisture in the spring

EDIT 05/12/16: CA and MX added, but only yield, production and harvested area data calculated
                The climate data require a subnational grid to aggregate onto, and so are skipped
                because the MIRCA data doesn't include them.
                
EDIT 06/02/16: fixed the cropland mask portion of the script so that variables could be binned
                in only the cropland regions.

EDIT 06/15/16: added the potential for soil water content (SWC) to be turned on. Intended for BR only

EDIT 07/05/16: added factors to convert everything to standard units of
    YIELD: kg/ha
    HARVESTED AREA: ha
    PRODUCTION: kg
    
EDIT: 11/11/16: fixed China conversion factor .;
"""

import numpy as np
import pandas as pd
import time
from scipy import signal
from scipy import ndimage
import netCDF4
import copy
start = time.clock()

def moving_average(a, smooth) :
    n=smooth*2+1
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

crops = ['maize']#
stORreg = 'state' #state or region
countries = ['MX']#
skipClim = 'T' #T/F
areas = ['all']#'all','AllArea'] # ,'halfPercMask' #all, AllArea or halfPercMask
GDD = 'GDD2'
dt = '' # '_dt' or ''
dep = 1 #for beginning of the flowering season soil moisture 0: 0-10cm, 1:10-40cm, 2:40-100, 3:100-200m
useSM = 'N'
notes = ''
prcpmn = 'all'
useSWC = 'N'

for area in areas:
    print(area)
    for crop in crops:
        print(crop)
        if crop=='springwheat':
            cropNam = 'w_sp'
        elif crop=='winterwheat':
            cropNam = 'w_win'
        else:
            cropNam = crop[0]
        
        #read in land-sea masks
        msk_gauss = netCDF4.Dataset('/Volumes/Data_Archive/Data/landSeaMask/gaussGrid/lsmask.nc')
        msk_1 = netCDF4.Dataset('/Volumes/Data_Archive/Data/landSeaMask/1.0Deg/lsmask.nc')
        LSmskGauss = msk_gauss.variables['land'][0,...]
        LSmsk1deg = msk_1.variables['mask'][0,...];
        
        #read in the data: spatial MIRCA country codes and FAOSTAT data
        cdasMircaCodes = np.load('/Volumes/Data_Archive/Data/CropCalendar/'+
                            'CCCmask/CDAS_grid/UnitCode/SubnatUnitCode.npy')
        gpccMircaCodes = np.load('/Volumes/Data_Archive/Data/CropCalendar/'+
                            'CCCmask/gpcc_grid/UnitCode/SubnatUnitCode.npy')
        noahMircaCodes = np.load('/Volumes/Data_Archive/Data/CropCalendar/'+
                            'CCCmask/Noah_grid/UnitCode/SubnatUnitCode.npy')
        sheffMircaCodes = np.load('/Volumes/Data_Archive/Data/CropCalendar/'+
                            'CCCmask/Sheffield_grid/UnitCode/SubnatUnitCode.npy')
                         
        #limit the analysis to major cropped area, all area or all cropped area
        if area =='all':
            CCC_gpcc = np.ones([gpccMircaCodes.shape[0],gpccMircaCodes.shape[1]])
            CCC_cdas = np.ones([cdasMircaCodes.shape[0],cdasMircaCodes.shape[1]])
        else:
            CCC_gpcc = np.load('/Volumes/Data_Archive/Results/globalCorrelatedRisks/simpleCropMaps/'+area+'/gpcc_grid/'+cropNam+'Map.npy')
            CCC_gpcc = np.ma.filled(CCC_gpcc,0)
            CCC_gpcc[CCC_gpcc==-99]=0 #set minor areas (code -99) equal to 0
            CCC_gpcc[CCC_gpcc>0]=1
            CCC_cdas = np.load('/Volumes/Data_Archive/Results/globalCorrelatedRisks/simpleCropMaps/'+area+'/CDAS_grid/'+cropNam+'Map.npy')
            CCC_cdas = np.ma.filled(CCC_cdas,0)
            CCC_cdas[CCC_cdas==-99]=0 #set minor areas (code -99) equal to 0
            CCC_cdas[CCC_cdas>0]=1
        
        cdasMircaCodes = cdasMircaCodes*CCC_cdas
        gpccMircaCodes = gpccMircaCodes*CCC_gpcc
        noahMircaCodes = noahMircaCodes*CCC_gpcc[:-30,:]
        sheffMircaCodes = sheffMircaCodes*CCC_gpcc
        
        #cdas corners need to be switched from 0 to 180 to match the nc grids
        cdasMircaCodes = np.append(cdasMircaCodes[:,97:],cdasMircaCodes[:,:97],axis=1)
        
        #limit the analysis to land areas, so that when you sum of subnational areas
        # you don't get values over the ocean                   
        cdasMircaCodes = cdasMircaCodes*LSmskGauss
        sheffMircaCodes = np.append(sheffMircaCodes[:,180:],sheffMircaCodes[:,:180],axis=1)  
        sheffMircaCodes = sheffMircaCodes*LSmsk1deg
        gpccMircaCodes = gpccMircaCodes*np.append(LSmsk1deg[:,180:],LSmsk1deg[:,:180],axis=1)  
        noahMircaCodes = noahMircaCodes*np.append(LSmsk1deg[:-30,180:],LSmsk1deg[:-30,:180],axis=1)  
        
        
        
        adminCodes = pd.read_table('/Volumes/Data_Archive/Data/CropCalendar/MIRCA2000/unit_name_state.txt',sep='\t')
        
        #re-parse the unit code names
        for n in range(0,np.size(adminCodes['UnitName'])):
            if np.size(adminCodes['UnitName'][n].split('_'))>1:
                newName = adminCodes['UnitName'][n].split('_')[1]
                adminCodes['UnitName'][n]=newName
        
        petNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/PET.nc')
        petNCenrm = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/PET_ensoRM.nc')
        vpdNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/VPD.nc')
        vpdNCenrm = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/VPD_ensoRM.nc')
        smNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/SoilM.nc')
        smNCenrm = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/SoilM_ensoRM.nc')
        gddNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/'+crop[0]+GDD+'.nc')
        gddNCenrm = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/'+crop[0]+GDD+'_ensoRM.nc')
        kddNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/'+crop[0]+'KDD.nc')
        kddNCenrm = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/'+crop[0]+'KDD_ensoRM.nc')
        sumpNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/sumP.nc')
        sumpNCenrm = netCDF4.Dataset('/Volumes/Data_Archive/Results/AmENSO/IntegratedVars'+dt+'/monthlyGlobal/sumP_ensoRM.nc')
        swc_ET0 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/referencepct.npy')
        swc_s1 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/south1pct.npy')
        swc_sS = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/south1pct.npy')
        swc_s2 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/south2pct.npy')
        swc_s3 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/south3pct.npy')
        swc_cer = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/cerradopct.npy')
        
        ENyrs = [1953, 1957, 1963, 1965, 1968, 1972, 1976, 1979, 1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009] #EN
        LNyrs = [1954, 1964, 1970, 1973, 1983, 1988, 1995, 1998, 2007, 2010] #LN0_2
        yearNC = gddNC.variables['year'][:]
        yearNCenrm = gddNCenrm.variables['year'][:]
        colNams = ['yield','yldAnom','yldAnomGau','yldAnomFD','yldAnomAbs','yldAnomSmooth5','yldAnomSmooth9',\
                    'yldAnomGauAbs','yldAnomSmooth5Abs','yldAnomSmooth9Abs','expectedYld','expectedYldSmooth5',\
                    'expectedYldSmooth9','expectedYldGau','Production','ProdAnomAbs',\
                    'Harvested_Area','GDD','KDD','KDD_sum','sumPan','sumPan_sum','sumP','sumP2','sumPanLag','PET','P_E',\
                    'VPD','GDD_ENrm','KDD_ENrm','sumPan_ENrm','sumP_ENrm','sumP2_ENrm','soilM',\
                    'soilMseed','sumPanLag_ENrm','soilM_ENrm','soilMseed_ENrm','PET_ENrm',\
                    'VPD_ENrm','EN','LN','LN2','SWC_ET0','SWC_s1','SWC_s2','SWC_s3','SWC_cer']
            
        mismatch = []
        
        for country in countries:
            stYr = [1950];endYr = [2010]
            yldFactor = 1.;prdFactor = 1.;HAfactor = 1.
            if 'sumpAnom_lag' in locals(): 
                del sumpAnom_lag
            #define the country parameters
            if country is 'US':
                dataFile = 'USDA_NASS'
                breakPt = 30 #number of years in to put the breakpoint for detrending
                
                if crop=='wheat':
                    prodFactor = 1000./36.7440 #Bu to tonnes to kg
                    yldFactor= (1000.*2.47)/36.7440 #Bu/Acre to Kg/ha
                    HAfactor = 1./2.47 #acres to ha
                    prcp_lag = ['Y'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([4,5,6])]#for precipitaiton, not python indexing
                    growMon = [np.array([4,5,6])] #3,7for temperature
                    smMon = [np.array([4])] #months for soil moisture
                    seedMonths = [np.array([-2,-1])] #october november seeding
                elif crop=='winterwheat':
                    prodFactor = 1000./36.7440
                    yldFactor = (1000.*2.47)/36.7440
                    HAfactor = 1./2.47
                    prcp_lag = ['Y'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([4,5,6])]#for precipitaiton, not python indexing
                    growMon = [np.array([4,5,6])] #3,7for temperature
                    smMon = [np.array([4])] #months for soil moisture
                    seedMonths = [np.array([-2,-1])] #october november seeding
                elif crop=='springwheat':
                    prodFactor = 1000./36.7440
                    yldFactor = (1000.*2.47)/36.7440
                    HAfactor = 1./2.47
                    prcp_lag = ['N'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([5,6,7])]#for precipitaiton, not python indexing
                    growMon = [np.array([5,6,7])] #3,7for temperature
                    smMon = [np.array([5])] #months for soil moisture
                    seedMonths = [np.array([-2,-1])] #october november seeding
                elif crop=='maize':
                    prodFactor = 1000./39.3680
                    yldFactor = (1000.*2.47)/39.3680
                    HAfactor = 1./2.47
                    prcp_lag = ['N'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([6,7,8])] #for precipitation
                    growMon = [np.array([6,7,8])] #9,10for temperature
                    smMon = [np.array([6])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([5])] #May seeding
                elif crop=='soy':
                    prodFactor = 1000./36.7440
                    yldFactor = (1000.*2.47)/36.7440
                    HAfactor = 1./2.47
                    prcp_lag = ['N'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([6,7,8])] #for precipitation
                    growMon = [np.array([6,7,8])] #9,10for temperature
                    smMon = [np.array([6])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([5])] #May seeding
                lag_months = [np.array([0,1,2,3])] #DJFM
                yldSt = [1950]; 
                yldEnd =[2013]
                
            elif country is 'AR':
                dataFile = 'AR_SIIA'
                prodFactor = 1000. #convert tonnes to kg
                yldFactor = 1000.#convert tonnes/ha to kg/ha
                breakPt = 30 #number of years in to put the breakpoint for detrending
                if crop=='soy':
                    prcp_lag = ['Y'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([13,14,15])] #for precipitation
                    growMon = [np.array([12,13,14])] #15,16 for temperature
                    smMon = [np.array([13])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([11,12])] 
                    yldSt = [1969]
                    yldEnd =[2013]
                elif crop=='maize':
                    prcp_lag = ['N'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([11,12,13])] 
                    growMon = [np.array([11,12,13])] #14,15
                    smMon = [np.array([11])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([9,10])] 
                    yldSt = [1969]
                    yldEnd =[2013]
                elif crop=='wheat':
                    prcp_lag = ['N'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([9,10,11])] 
                    growMon = [np.array([9,10,11])] #8,12
                    smMon = [np.array([9])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([6,7])] #June/Jul seeding
                    yldSt = [1969]
                    yldEnd =[2013]
                lag_months = [np.array([9,10,11,12])] #SOND
        
            elif country is 'BR':
                useSWC = 'Y'
                dataFile = 'BR_CONAB'
                prodFactor = 1000.*1000 #convert thousands of tonnes to kg
                yldFactor = 1.#strangely, BR claims to be in (Em kg/ha) but is clearly (kg/ha)
                HAfactor = 1000.#convert thousands of ha to ha
                lag_months = [np.array([9,10,11,12])] #SOND  
                if crop=='soy':
                    prcp_lag = ['Y'] #whether to include lagged precip in the regression
                    swc_cer_mo = [np.array([-4,-3,-2,-1,0])]#SONDJ
                    swc_s_mo = [np.array([-3,-2,-1,0,1,2])]#ONDJFM
                    prcpMon = [np.array([-3,-2,-1,0,1,2])] #[1,2,3]
                    growMon = [np.array([-3,-2,-1,0,1,2])] #
                    smMon = [np.array([1])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([-1,0])] 
                    yldSt = [1978]
                    yldEnd = [2013]
                    breakPt = 7 #number of years in to put the breakpoint for detrending
                elif (crop=='maize1')|(crop=='maize'):
                    prcp_lag = ['N'] #whether to include lagged precip in the regression
                    swc_cer_mo = [np.array([13,14,15,16])]#[np.array([14,15])]#
                    swc_s_mo = [np.array([10,11,12,13,14,15])]#ONDJFM
                    prcpMon = [np.array([9,10,11,12,13,14,15])] #[11,12,13]
                    growMon = [np.array([9,10,11,12,13,14,15])] #
                    smMon = [np.array([11])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([9])] 
                    yldSt = [1978]
                    yldEnd = [2013]
                    breakPt = 12
                elif crop=='maize2':
                    prcp_lag = ['Y'] #whether to include lagged precip in the regression
                    swc_cer_mo = [np.array([13,14,15,16])]#[np.array([14,15])]#
                    swc_s_mo = [np.array([10,11,12,13,14,15])]#ONDJFM
                    prcpMon = [np.array([13,14,15,16])] #14,15,16
                    growMon = [np.array([13,14,15,16])] #
                    smMon = [np.array([14])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([13])] 
                    yldSt =  [1978]
                    yldEnd = [2013]
                    breakPt = 0
                elif crop=='wheat':
                    useSWC = 'N'
                    prcp_lag = ['N'] #whether to include lagged precip in the regression
                    prcpMon = [np.array([9,10,11])] 
                    growMon = [np.array([9,10,11])]#8,12
                    smMon = [np.array([9])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([5,6])] #June seeding
                    yldSt = [1978]
                    yldEnd = [2013]
                    breakPt = 7
                    
            elif country is 'MX':
                prodFactor = 1000. #convert tonnes to kg:
                yldFactor = 1000. #convert tonnes/ha to kg/ha             
                prcp_lag = ['N'] #whether to include lagged precip in the regression
                dataFile = 'MX_SAGARPA'
                yldSt = [1980]
                yldEnd =[2012]
                lag_months = [np.array([])]                    
                prcpMon = [np.array([5,6,7])] 
                growMon = [np.array([5,6,7])]
                smMon = [np.array([5])] #beginning of flowering months for soil moisture
                seedMonths = [np.array([-1,0])] #ND planting?
                breakPt = 0 #number of years in to put the breakpoint for detrending   
               
            elif country is 'CA':
                prcp_lag = ['N'] #whether to include lagged precip in the regression
                dataFile = 'CA_CANSIM'
                prodFactor = 1000. #convert tonnes to kg
                yldFactor = 1000. #convert tonnes/ha to kg/ha   
                yldSt = [1950]
                yldEnd =[2013]
                lag_months = [np.array([2,3,4])]                    
                prcpMon = [np.array([5,6,7])] 
                growMon = [np.array([5,6,7])]
                smMon = [np.array([5])] #beginning of flowering months for soil moisture
                seedMonths = [np.array([4,5])] #June seeding
                breakPt = 30 #number of years in to put the breakpoint for detrending   
            
            elif country is 'AU':
                prcp_lag = ['N'] #whether to include lagged precip in the regression
                dataFile = 'AUS_ABS'
                prodFactor = 1000. #convert tonnes to kg:
                yldFactor = 1000. #convert tonnes/ha to kg/ha      
                yldSt = [1950]
                yldEnd =[2011]
                lag_months = [np.array([6,7,8])]                    
                prcpMon = [np.array([9,10,11])] 
                growMon = [np.array([9,10,11])]
                smMon = [np.array([9])] #beginning of flowering months for soil moisture
                seedMonths = [np.array([8,9])]
                breakPt = 30 #number of years in to put the breakpoint for detrending   

            elif country is 'CHN':
                prcp_lag = ['N'] #whether to include lagged precip in the regression
                dataFile = 'CHN'
                yldSt = [1949]
                yldEnd =[2013]
                prodFactor = 10000*1000. #10k tonnes to kg
                HAfactor = 666.7 #10k "management units" (mu) to ha
                yldFactor = (10000/666.66) # kg/mu to kg/ha
                if (crop=='wheat')|(crop=='maize')|(crop=='winterwheat'):
                    lag_months = [np.array([1,2,3])]                    
                    prcpMon = [np.array([4,5,6])] 
                    growMon = [np.array([4,5,6])]
                    smMon = [np.array([4])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([3,4])] 
                elif crop=='soy':
                    lag_months = [np.array([4,5,6])]                    
                    prcpMon = [np.array([7,8,9])] 
                    growMon = [np.array([7,8,9])]
                    smMon = [np.array([7])] #beginning of flowering months for soil moisture
                    seedMonths = [np.array([6,7])] 
                breakPt = 30 #number of years in to put the breakpoint for detrending   
                        
            
                        
            
            endYr[0]=endYr[0]-1    
            yrNCind = (yearNC<=endYr)&(yearNC>=stYr)
            yrNCindenrm = (yearNCenrm<=endYr)&(yearNCenrm>=stYr)
            
            #only deal with the climate data if need be
            if skipClim !='T':
            #calculate the SWC
            ########################## SWC ##########################
                if useSWC is 'Y':
                    ind_s = np.repeat((stYr[0] - 1950+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(swc_s_mo))+np.tile(swc_s_mo[0]-1,endYr[0]+1-stYr[0])
                    ind_cer = np.repeat((stYr[0] - 1950+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(swc_cer_mo))+np.tile(swc_cer_mo[0]-1,endYr[0]+1-stYr[0])
                    swc_s1 = swc_s1[ind_s,...];swc_s2 = swc_s2[ind_s,...];swc_s3 = swc_s3[ind_s,...];
                    swc_cer = swc_cer[ind_cer,...];swc_ET0 = swc_ET0[ind_s,...]
                    swc_s1=np.reshape(swc_s1,[np.size(swc_s_mo),endYr[0]+1-stYr[0],swc_s1.shape[1],swc_s1.shape[2]],order='F')
                    swc_s2=np.reshape(swc_s2,[np.size(swc_s_mo),endYr[0]+1-stYr[0],swc_s2.shape[1],swc_s2.shape[2]],order='F')
                    swc_s3=np.reshape(swc_s3,[np.size(swc_s_mo),endYr[0]+1-stYr[0],swc_s3.shape[1],swc_s3.shape[2]],order='F')
                    swc_cer=np.reshape(swc_cer,[np.size(swc_cer_mo),endYr[0]+1-stYr[0],swc_cer.shape[1],swc_cer.shape[2]],order='F')
                    swc_ET0=np.reshape(swc_ET0,[np.size(swc_s_mo),endYr[0]+1-stYr[0],swc_ET0.shape[1],swc_ET0.shape[2]],order='F')
                    swc_s1 = np.nanmean(swc_s1,axis=0)
                    swc_s2 = np.nanmean(swc_s2,axis=0)
                    swc_s3 = np.nanmean(swc_s3,axis=0)
                    swc_cer = np.nanmean(swc_cer,axis=0)
                    swc_ET0 = np.nanmean(swc_ET0,axis=0)
                    
            #calculate the GDD, KDD, sumP, and soil moisture
            ##################### soil moisture #####################
                if useSM is 'Y':
                    ind = np.repeat((stYr[0] - smNC.variables['year'][0]+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(smMon))+np.tile(smMon[0]-1,endYr[0]+1-stYr[0])
                    sm = smNC.variables['SoilM'][:,:,dep,:,:] # pull the soil moisture for the specified depth
                    sm = np.swapaxes(sm,0,1)        
                    sm = np.reshape(sm,[sm.shape[0]*sm.shape[1],sm.shape[2],sm.shape[3]],order='F')
                    sm = sm[ind,...]; sm=np.reshape(sm,[np.size(smMon),endYr[0]+1-stYr[0],sm.shape[1],sm.shape[2]],order='F')
                    sm = np.ma.filled(sm); sm = sm[:,:,::-1,...];
                    sm = np.nanmean(sm,axis=0)
                    smAnom = sm - np.nanmean(sm,axis=0) 
                    
                    smENrm = smNCenrm.variables['SoilM'][:,:,dep,:,:] # pull the soil moisture for the specified depth
                    smENrm = np.swapaxes(smENrm,0,1)        
                    smENrm = np.reshape(smENrm,[smENrm.shape[0]*smENrm.shape[1],smENrm.shape[2],smENrm.shape[3]],order='F')
                    smENrm = smENrm[ind,...]; smENrm=np.reshape(smENrm,[np.size(smMon),endYr[0]+1-stYr[0],smENrm.shape[1],smENrm.shape[2]],order='F')
                    smENrm = np.ma.filled(smENrm); smENrm = smENrm[:,:,::-1,...]
                    smENrm = np.nansum(smENrm,axis=0)
                    smAnomENrm = smENrm - np.nanmean(sm,axis=0)
            
                    #do the same for the planting soil moisture
                    ind = np.repeat((stYr[0] - smNC.variables['year'][0]+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(seedMonths))+np.tile(seedMonths[0],endYr[0]+1-stYr[0])
                    smSeed = smNC.variables['SoilM'][:,:,dep,:,:] # pull the soil moisture for the specified depth
                    smSeed = np.swapaxes(smSeed,0,1)        
                    smSeed = np.reshape(smSeed,[smSeed.shape[0]*smSeed.shape[1],smSeed.shape[2],smSeed.shape[3]],order='F')
                    smSeed = smSeed[ind,...]; smSeed=np.reshape(smSeed,[np.size(seedMonths),endYr[0]+1-stYr[0],smSeed.shape[1],smSeed.shape[2]],order='F')
                    smSeed = np.ma.filled(smSeed); smSeed = smSeed[:,:,::-1,...]; 
                    smSeed = np.nanmean(smSeed,axis=0)
                    smSeedAnom = smSeed - np.nanmean(smSeed,axis=0) 
                    
                    smSeedENrm = smNCenrm.variables['SoilM'][:,:,dep,:,:] # pull the soil moisture for the specified depth
                    smSeedENrm = np.swapaxes(smSeedENrm,0,1)        
                    smSeedENrm = np.reshape(smSeedENrm,[smSeedENrm.shape[0]*smSeedENrm.shape[1],smSeedENrm.shape[2],smSeedENrm.shape[3]],order='F')
                    smSeedENrm = smSeedENrm[ind,...]; smSeedENrm=np.reshape(smSeedENrm,[np.size(seedMonths),endYr[0]+1-stYr[0],smSeedENrm.shape[1],smSeedENrm.shape[2]],order='F')
                    smSeedENrm = np.ma.filled(smSeedENrm); smSeedENrm = smSeedENrm[:,:,::-1,...]
                    smSeedENrm = np.nansum(smSeedENrm,axis=0)
                    smSeedAnomENrm = smSeedENrm - np.nanmean(smSeed,axis=0)     
                    
                ##################### PET #####################
                ind = np.repeat((stYr[0] - petNC.variables['year'][0]+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(growMon))+np.tile(growMon[0]-1,endYr[0]+1-stYr[0])
                pet = petNC.variables['PET'][:] # pull the soil moisture for the specified depth
                pet = np.swapaxes(pet,0,1)        
                pet = np.reshape(pet,[pet.shape[0]*pet.shape[1],pet.shape[2],pet.shape[3]],order='F')
                pet = pet[ind,...]; pet=np.reshape(pet,[np.size(growMon),endYr[0]+1-stYr[0],pet.shape[1],pet.shape[2]],order='F')
                pet = np.ma.filled(pet); pet = pet[:,:,::-1,...];
                sumPET = np.nansum(pet,axis=0) #calculate the sum for P-E later
                pet = np.nanmean(pet,axis=0); pet=pet[:,:-30,:]
                petAnom = pet - np.nanmean(pet,axis=0) 
                
                petENrm = petNCenrm.variables['PET'][:] # pull the soil moisture for the specified depth
                petENrm = np.swapaxes(petENrm,0,1)        
                petENrm = np.reshape(petENrm,[petENrm.shape[0]*petENrm.shape[1],petENrm.shape[2],petENrm.shape[3]],order='F')
                petENrm = petENrm[ind,...]; petENrm=np.reshape(petENrm,[np.size(growMon),endYr[0]+1-stYr[0],petENrm.shape[1],petENrm.shape[2]],order='F')
                petENrm = np.ma.filled(petENrm); petENrm = petENrm[:,:,::-1,...]; 
                petENrm = np.nanmean(petENrm,axis=0)
                petAnomENrm = petENrm - np.nanmean(pet,axis=0)
                ##################### VPD #####################
                vpd = vpdNC.variables['VPD'][:] # pull the soil moisture for the specified depth
                vpd = np.swapaxes(vpd,0,1)        
                vpd = np.reshape(vpd,[vpd.shape[0]*vpd.shape[1],vpd.shape[2],vpd.shape[3]],order='F')
                vpd = vpd[ind,...]; vpd=np.reshape(vpd,[np.size(growMon),endYr[0]+1-stYr[0],vpd.shape[1],vpd.shape[2]],order='F')
                vpd = np.ma.filled(vpd); vpd = vpd[:,:,::-1,...]; 
                vpd = np.nanmean(vpd,axis=0); vpd=vpd[:,:-30,:]
                vpdAnom = vpd - np.nanmean(vpd,axis=0) 
                
                vpdENrm = vpdNCenrm.variables['VPD'][:] # pull the soil moisture for the specified depth
                vpdENrm = np.swapaxes(vpdENrm,0,1)        
                vpdENrm = np.reshape(vpdENrm,[vpdENrm.shape[0]*vpdENrm.shape[1],vpdENrm.shape[2],vpdENrm.shape[3]],order='F')
                vpdENrm = vpdENrm[ind,...]; vpdENrm=np.reshape(vpdENrm,[np.size(growMon),endYr[0]+1-stYr[0],vpdENrm.shape[1],vpdENrm.shape[2]],order='F')
                vpdENrm = np.ma.filled(vpdENrm); vpdENrm = vpdENrm[:,:,::-1,...]; 
                vpdENrm = np.nanmean(vpdENrm,axis=0)
                vpdAnomENrm = vpdENrm - np.nanmean(vpd,axis=0)
                ##################### GDD #####################
                ind = np.repeat((stYr[0] - gddNC.variables['year'][0]+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(growMon))+np.tile(growMon[0]-1,endYr[0]+1-stYr[0])    
                gdd = gddNC.variables[GDD][:] # pull the soil moisture for the specified depth
                gdd = np.swapaxes(gdd,0,1)        
                gdd = np.reshape(gdd,[gdd.shape[0]*gdd.shape[1],gdd.shape[2],gdd.shape[3]],order='F')
                gdd = gdd[ind,...]; gdd=np.reshape(gdd,[np.size(growMon),endYr[0]+1-stYr[0],gdd.shape[1],gdd.shape[2]],order='F')
                gdd = np.ma.filled(gdd); 
                gdd = np.nansum(gdd,axis=0)
                gddAnom = gdd - np.nanmean(gdd,axis=0) 
                
                indENrm = np.repeat((stYr[0] - gddNCenrm.variables['year'][0]+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(growMon))+np.tile(growMon[0]-1,endYr[0]+1-stYr[0])
                gddENrm = gddNCenrm.variables[GDD][:] # pull the soil moisture for the specified depth
                gddENrm = np.swapaxes(gddENrm,0,1)        
                gddENrm = np.reshape(gddENrm,[gddENrm.shape[0]*gddENrm.shape[1],gddENrm.shape[2],gddENrm.shape[3]],order='F')
                gddENrm = gddENrm[indENrm,...]; gddENrm=np.reshape(gddENrm,[np.size(growMon),endYr[0]+1-stYr[0],gddENrm.shape[1],gddENrm.shape[2]],order='F')
                gddENrm = np.ma.filled(gddENrm);
                gddENrm = np.nansum(gddENrm,axis=0)
                gddAnomENrm = gddENrm - np.nanmean(gdd,axis=0)
                ##################### KDD #####################
                kdd = kddNC.variables['KDD'][:]
                kdd = np.swapaxes(kdd,0,1)        
                kdd = np.reshape(kdd,[kdd.shape[0]*kdd.shape[1],kdd.shape[2],kdd.shape[3]],order='F')
                kdd = kdd[ind,...]; kdd=np.reshape(kdd,[np.size(growMon),endYr[0]+1-stYr[0],kdd.shape[1],kdd.shape[2]],order='F')
                kdd = np.ma.filled(kdd); 
                kdd = np.nansum(kdd,axis=0)
                kddAnom = kdd - np.nanmean(kdd,axis=0) 
                
                kddENrm = kddNCenrm.variables['KDD'][:] # pull the soil moisture for the specified depth
                kddENrm = np.swapaxes(kddENrm,0,1)        
                kddENrm = np.reshape(kddENrm,[kddENrm.shape[0]*kddENrm.shape[1],kddENrm.shape[2],kddENrm.shape[3]],order='F')
                kddENrm = kddENrm[indENrm,...]; kddENrm=np.reshape(kddENrm,[np.size(growMon),endYr[0]+1-stYr[0],kddENrm.shape[1],kddENrm.shape[2]],order='F')
                kddENrm = np.ma.filled(kddENrm); 
                kddENrm = np.nansum(kddENrm,axis=0)
                kddAnomENrm = kddENrm - np.nanmean(kdd,axis=0)
                ##################### sumP #####################
                ind = np.repeat((stYr[0] - sumpNC.variables['year'][0]+np.array(range(endYr[0]+1-stYr[0])))*12,np.size(prcpMon))+np.tile(prcpMon[0]-1,endYr[0]+1-stYr[0])
                sump = sumpNC.variables['sumP'][:] # pull the soil moisture for the specified depth
                sump = np.swapaxes(sump,0,1)        
                sump = np.reshape(sump,[sump.shape[0]*sump.shape[1],sump.shape[2],sump.shape[3]],order='F')
                sump = sump[ind,...]; sump=np.reshape(sump,[np.size(prcpMon),endYr[0]+1-stYr[0],sump.shape[1],sump.shape[2]],order='F')
                sump = np.ma.filled(sump); 
                sump = np.nansum(sump,axis=0)
                P_E = sump-sumPET #calculate P - E for the growing season
                sumpAnom = sump - np.nanmean(sump,axis=0) 
                
                sumpENrm = sumpNCenrm.variables['sumP'][:] # pull the soil moisture for the specified depth
                sumpENrm = np.swapaxes(sumpENrm,0,1)        
                sumpENrm = np.reshape(sumpENrm,[sumpENrm.shape[0]*sumpENrm.shape[1],sumpENrm.shape[2],sumpENrm.shape[3]],order='F')
                sumpENrm = sumpENrm[ind,...]; sumpENrm=np.reshape(sumpENrm,[np.size(prcpMon),endYr[0]+1-stYr[0],sumpENrm.shape[1],sumpENrm.shape[2]],order='F')
                sumpENrm = np.ma.filled(sumpENrm); 
                sumpENrm = np.nansum(sumpENrm,axis=0)
                sumpAnomENrm = sumpENrm - np.nanmean(sump,axis=0)               
                if (prcp_lag[0]=='Y'):
                    ind = np.repeat((stYr[0] - sumpNC.variables['year'][0]+np.array(range(endYr[0]+1-stYr[0])))*12,4)+np.tile(lag_months[0]-1,endYr[0]+1-stYr[0])
                    sumpLag = sumpNC.variables['sumP'][:] # pull the soil moisture for the specified depth
                    sumpLag = np.swapaxes(sumpLag,0,1)        
                    sumpLag = np.reshape(sumpLag,[sumpLag.shape[0]*sumpLag.shape[1],sumpLag.shape[2],sumpLag.shape[3]],order='F')
                    sumpLag = sumpLag[ind,...]; sumpLag=np.reshape(sumpLag,[4,endYr[0]+1-stYr[0],sumpLag.shape[1],sumpLag.shape[2]],order='F')
                    sumpLag = np.ma.filled(sumpLag); 
                    sumpLag = np.nansum(sumpLag,axis=0)
                    sumpAnom_lag = sumpLag - np.nanmean(sumpLag,axis=0) 
                    
                    sumpLagENrm = sumpNCenrm.variables['sumP'][:] # pull the soil moisture for the specified depth
                    sumpLagENrm = np.swapaxes(sumpLagENrm,0,1)        
                    sumpLagENrm = np.reshape(sumpLagENrm,[sumpLagENrm.shape[0]*sumpLagENrm.shape[1],sumpLagENrm.shape[2],sumpLagENrm.shape[3]],order='F')
                    sumpLagENrm = sumpLagENrm[ind,...]; sumpLagENrm=np.reshape(sumpLagENrm,[4,endYr[0]+1-stYr[0],sumpLagENrm.shape[1],sumpLagENrm.shape[2]],order='F')
                    sumpLagENrm = np.ma.filled(sumpLagENrm); 
                    sumpLagENrm = np.nansum(sumpLagENrm,axis=0)
                    sumpAnom_lagENrm = sumpLagENrm - np.nanmean(sumpLag,axis=0)
            
                    #fill the first year (which has no lag) with nans
                    #sumpAnom_lagENrm = np.append(np.zeros([1,sumpAnom.shape[1],sumpAnom.shape[2]])*np.nan,sumpAnom_lagENrm,axis=0)        
                    #if area == 'major':
                    #    sumpAnom_lagENrm = sumpAnom_lagENrm*CCC_gpcc
        
            #read in the yield data
            CNTRYyld = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/Yield/HistoricalData/'+dataFile+'/'+country+'_'+crop+'_'+stORreg+'.csv',
                                          header=0); CNTRYyld = CNTRYyld.sort()
            #read in the harvested area data
            CNTRYha = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/HarvestedArea/HistoricalData/'+dataFile+'/'+country+'_'+crop+'_'+stORreg+'.csv',
                                          header=0); CNTRYha = CNTRYha.sort()
            #read in the total production data
            CNTRYp = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/production/HistoricalData/'+dataFile+'/'+country+'_'+crop+'_'+stORreg+'.csv',
                                          header=0)
            CNTRYp = CNTRYp.sort();
            #NaNs coded as 0 in the data
            CNTRYyld[CNTRYyld==0]=np.nan
            CNTRYha[CNTRYha==0]=np.nan
            
            states = CNTRYyld.columns
        
            #create the pandas dataframe to fill, use a multi index of [year, state]
            # To get all state/year combinations 
            index = pd.MultiIndex.from_product([yearNC,states], names=['year', 'state'])
            df = pd.DataFrame(index=index, columns=colNams)    
            
            #linearly detrend each state yields, then 
            for stNam in states:
                stYield = CNTRYyld[str(yldSt[0]):str(yldEnd[0])][stNam]
                HA = np.array(CNTRYha[str(yldSt[0]):str(yldEnd[0])][stNam])
                Prod = np.array(CNTRYp[str(yldSt[0]):str(yldEnd[0])][stNam])
                Prod=Prod[Prod!=0]; Prod=Prod[~np.isnan(Prod)]; HA = HA[~np.isnan(HA)]; stYield = stYield.dropna();
                stYield = stYield*yldFactor
                Prod = Prod*prodFactor
                HA = HA*HAfactor 
                if (len(stYield) <= breakPt)|(len(Prod) <= breakPt)|(len(stYield)<9):
                    continue
                if np.nanmax(100*(stYield.values*HA - Prod)/Prod)>2.:
                    print(stNam)
                    print('something is wrong, check units of yield, HA and production')
                    print(np.nanmax(100*(stYield.values*HA - Prod)/Prod))
                        #break
                        #if the difference in yield*HA and production is more than 2%, then something is wrong
                
                stYrs = stYield.index.year
                stProdAbs = np.array(signal.detrend(Prod,bp=breakPt))
                stYldAbs = np.array(signal.detrend(stYield,bp=breakPt))
                stYld = np.array(signal.detrend(stYield,bp=breakPt)/
                        (stYield-signal.detrend(stYield,bp=breakPt))) #%yield anom = (yld obs - yld exp)/yld exp
                expYld = (stYield-signal.detrend(stYield,bp=breakPt))
                stYldFD = np.array(stYield[1:])-np.array(stYield[:-1]) 
                stYldSm9 = np.array(stYield[4:-4]-moving_average(np.array(stYield),4))/moving_average(np.array(stYield),4)
                stYldSm5 = np.array(stYield[2:-2]-moving_average(np.array(stYield),2))/moving_average(np.array(stYield),2)   
                stYldSm9Abs= np.array(stYield[4:-4]-moving_average(np.array(stYield),4))
                stYldSm5Abs = np.array(stYield[2:-2]-moving_average(np.array(stYield),2))
                stYldGau = (stYield.values-ndimage.filters.gaussian_filter1d(stYield.values,3))/ndimage.filters.gaussian_filter1d(stYield.values,3)
                select = list(zip(stYrs,np.repeat(stNam,stYrs.size)))
                df.loc[select,'yield']=np.array(stYield)
                df.loc[select,'Production']=Prod
                df.loc[select,'ProdAnomAbs']=stProdAbs
                df.loc[select,'yldAnom']=stYld
                df.loc[select,'yldAnomGau']=stYldGau
                df.loc[select,'expectedYld']=np.array(expYld)
                df.loc[select[2:-2],'expectedYldSmooth5']=moving_average(np.array(stYield),2)
                df.loc[select[4:-4],'expectedYldSmooth9']=moving_average(np.array(stYield),4)
                df.loc[select,'expectedYldGau']=ndimage.filters.gaussian_filter1d(stYield.values,3)
                df.loc[select,'yldAnomAbs']=stYldAbs
                df.loc[select,'yldAnomGauAbs']=(stYield.values-ndimage.filters.gaussian_filter1d(stYield.values,3))
                df.loc[select[2:-2],'yldAnomSmooth5']=stYldSm5
                df.loc[select[4:-4],'yldAnomSmooth9']=stYldSm9
                df.loc[select[2:-2],'yldAnomSmooth5Abs']=stYldSm5Abs
                df.loc[select[4:-4],'yldAnomSmooth9Abs']=stYldSm9Abs
                if np.shape(select)[0]>1:
                    df.loc[select[1:],'yldAnomFD']=stYldFD
                df.loc[select,'Harvested_Area']=HA            
        
                ENselect = list(zip(np.array(ENyrs),np.repeat(stNam,np.size(ENyrs))))
                ENselect_1 = list(zip(np.array(ENyrs)-1,np.repeat(stNam,np.size(ENyrs))))
                ENselect1 = list(zip(np.array(ENyrs)+1,np.repeat(stNam,np.size(ENyrs))))
                LNselect = list(zip(np.array(LNyrs),np.repeat(stNam,np.size(LNyrs))))
                LNselect_1 = list(zip(np.array(LNyrs)-1,np.repeat(stNam,np.size(LNyrs))))
                LNselect1 = list(zip(np.array(LNyrs)+1,np.repeat(stNam,np.size(LNyrs))))
                LNselect2 = list(zip(np.array(LNyrs)+2,np.repeat(stNam,np.size(LNyrs))))
                df.loc[ENselect,'EN']=0;df.loc[ENselect_1,'EN']=-1;df.loc[ENselect1,'EN']=1
                df.loc[LNselect,'LN']=0;df.loc[LNselect_1,'LN']=-1;
                df.loc[LNselect,'LN2']=0;df.loc[LNselect_1,'LN2']=-1;#repeat for when LN2 == LN-1 (so it is overwritten)
                if yearNC[yrNCind][-1] == 2010: #necessary because no 2011 soilm data
                    df.loc[LNselect1[:-1],'LN']=1
                    df.loc[LNselect1[:-1],'LN2']=1
                    df.loc[LNselect2[:-1],'LN2']=2
                else:
                    df.loc[LNselect1,'LN']=1
                    df.loc[LNselect1,'LN2']=1
                    df.loc[LNselect2,'LN2']=2
                    
                #Skip the rest of the climate analyses if indicated        
                if skipClim=='T':continue        
                
                noahMircaCodesTemp = copy.deepcopy(noahMircaCodes)
                cdasMircaCodesTemp = copy.deepcopy(cdasMircaCodes)
                sheffMircaCodesTemp = copy.deepcopy(sheffMircaCodes[:-30,:])
                gpccMircaCodesTemp = copy.deepcopy(gpccMircaCodes)
                code = adminCodes['UnitCode'][adminCodes['UnitName']==stNam.title()]
                if code.size==0:
                    mismatch.append(stNam)
                    continue
                for num in code:
                    noahMircaCodesTemp[noahMircaCodes==num]=1
                    cdasMircaCodesTemp[cdasMircaCodes==num]=1
                    sheffMircaCodesTemp[sheffMircaCodes[:-30,:]==num]=1
                    gpccMircaCodesTemp[gpccMircaCodes==num]=1
                noahMircaCodesTemp[noahMircaCodesTemp!=1]=0
                cdasMircaCodesTemp[cdasMircaCodesTemp!=1]=0
                sheffMircaCodesTemp[sheffMircaCodesTemp!=1]=0
                gpccMircaCodesTemp[gpccMircaCodesTemp!=1]=0     
                
                if np.array(CNTRYyld[str(yldSt[0]):str(yldEnd[0])][stNam]).dtype == 'O':
                    continue #this skips columns with missing values
                #create the index to use in selecting 
                select = list(zip(yearNC[yrNCind],np.repeat(stNam,yearNC[yrNCind].size)))
                
                if useSWC is 'Y':
                    swcET0AnomTS = np.nanmean(swc_ET0[:,gpccMircaCodesTemp==1],axis=1)
                    df.loc[select,'SWC_ET0']=swcET0AnomTS
                    swcS1AnomTS = np.nanmean(swc_s1[:,gpccMircaCodesTemp==1],axis=1)
                    df.loc[select,'SWC_s1']=swcS1AnomTS
                    swcS2AnomTS = np.nanmean(swc_s2[:,gpccMircaCodesTemp==1],axis=1)
                    df.loc[select,'SWC_s2']=swcS2AnomTS
                    swcS3AnomTS = np.nanmean(swc_s3[:,gpccMircaCodesTemp==1],axis=1)
                    df.loc[select,'SWC_s3']=swcS3AnomTS 
                    swcCerAnomTS = np.nanmean(swc_cer[:,gpccMircaCodesTemp==1],axis=1)
                    df.loc[select,'SWC_cer']=swcCerAnomTS 
                if useSM is 'Y':
                    smAnomTS = np.nanmean(smAnom[:,noahMircaCodesTemp==1],axis=1)
                    df.loc[select,'soilM']=smAnomTS
                    smAnomTS_ENrm = np.nanmean(smAnomENrm[:,noahMircaCodesTemp==1],axis=1)
                    df.loc[select,'soilM_ENrm']=smAnomTS_ENrm
                    smSeedAnomTS_ENrm = np.nanmean(smSeedAnomENrm[:,noahMircaCodesTemp==1],axis=1)
                    df.loc[select,'soilMseed_ENrm']=smSeedAnomTS_ENrm
                    smSeedAnomTS = np.nanmean(smSeedAnom[:,noahMircaCodesTemp==1],axis=1)
                    df.loc[select,'soilMseed']=smSeedAnomTS
                gddAnomTS = np.nanmean(gddAnom[:,cdasMircaCodesTemp==1],axis=1)
                df.loc[select,'GDD']=gddAnomTS
                petAnomTS = np.nanmean(petAnom[:,sheffMircaCodesTemp==1],axis=1)
                df.loc[select,'PET']=petAnomTS
                vpdAnomTS = np.nanmean(vpdAnom[:,sheffMircaCodesTemp==1],axis=1)
                df.loc[select,'VPD']=vpdAnomTS
                kddAnomTS = np.nanmean(kddAnom[:,cdasMircaCodesTemp==1],axis=1)
                kddAnomTS_sum = np.nansum(kddAnom[:,cdasMircaCodesTemp==1],axis=1)
                df.loc[select,'KDD']=kddAnomTS
                df.loc[select,'KDD_sum']=kddAnomTS_sum
                sumpAnomTS = np.nanmean(sumpAnom[:,gpccMircaCodesTemp==1],axis=1)
                sumpAnomTS_sum = np.nansum(sumpAnom[:,gpccMircaCodesTemp==1],axis=1)
                df.loc[select,'sumPan']=sumpAnomTS
                df.loc[select,'sumPan_sum']=sumpAnomTS_sum
                P_E_TS = np.nanmean(P_E[:,gpccMircaCodesTemp==1],axis=1)
                df.loc[select,'P_E']=P_E_TS
                sumpTS = np.nanmean(sump[:,gpccMircaCodesTemp==1],axis=1)
                df.loc[select,'sumP']=sumpTS
                df.loc[select,'sumP2']=sumpAnomTS*sumpAnomTS
        
                gddAnomTSENrm = np.nanmean(gddAnomENrm[:,cdasMircaCodesTemp==1],axis=1)
                df.loc[select,'GDD_ENrm']=gddAnomTSENrm
                petAnomTSENrm = np.nanmean(petAnomENrm[:,sheffMircaCodesTemp==1],axis=1)
                df.loc[select,'PET_ENrm']=petAnomTSENrm
                vpdAnomTSENrm = np.nanmean(vpdAnomENrm[:,sheffMircaCodesTemp==1],axis=1)
                df.loc[select,'VPD_ENrm']=vpdAnomTSENrm
                kddAnomTSENrm = np.nanmean(kddAnomENrm[:,cdasMircaCodesTemp==1],axis=1)
                df.loc[select,'KDD_ENrm']=kddAnomTSENrm
                sumpAnomTSENrm = np.nanmean(sumpAnomENrm[:,gpccMircaCodesTemp==1],axis=1)
                df.loc[select,'sumPan_ENrm']=sumpAnomTSENrm
                sumpTSENrm = np.nanmean(sumpENrm[:,gpccMircaCodesTemp==1],axis=1)
                df.loc[select,'sumP_ENrm']=sumpTSENrm
                df.loc[select,'sumP2_ENrm']=sumpAnomTSENrm*sumpAnomTSENrm
                if prcp_lag[0]=='Y':
                    sumpAnom_lagTS = np.nanmean(sumpAnom_lag[:,gpccMircaCodesTemp==1],axis=1)
                    df.loc[select,'sumPanLag']=sumpAnom_lagTS
                    sumpAnom_lagTSENrm = np.nanmean(sumpAnom_lagENrm[:,gpccMircaCodesTemp==1],axis=1)
                    df.loc[select,'sumPanLag_ENrm']=sumpAnom_lagTSENrm
        
                
            df.to_csv('/Volumes/Data_Archive/Results/AmENSO_yld/AmENSO_'+country+'_'+crop+'_'+area+GDD+dt+prcpmn+stORreg+notes+'.csv',',')
                
elapsed = time.clock()-start
print(elapsed)