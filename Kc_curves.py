# -*- coding: utf-8 -*-
"""
Take the soil water balance objects created by "soilWaterBalance.py" and 
correlate them to crop yield data

Created on Mon May  9 15:42:46 2016

@author: westonanderson
"""
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import copy
import netCDF4
from mpl_toolkits.basemap import Basemap
#plt.ioff()

#define the season to use
cropNam = 'maize' #soySouth, soyCerr, maize1, maize2, maize
lag = '' #lag or ''
SMcorr = 'T' # T/F
smData = 'IGBP' #Dunne or IGBP
yldAnom = 'yldAnom'
area = 'all' #allCrop or 'major' or 'all' or 'MODIScrop'
SMperc = 'T' #T/F
lonMin = -57; lonMax = -47; latMin = -15; latMax = -25 
years = [1950,2011];numYrs = years[1]-years[0]
ENyrs = np.array([1953, 1957, 1963, 1965, 1968, 1972, 1976, 1979, 1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009]) #EN
#LNyrs = np.array([1954, 1964, 1970, 1973, 1983, 1988, 1995, 1998, 2007, 2010]) #LN0 yrs
LNyrs = np.array([1983, 1984, 1988, 1995, 1996, 1998, 1999, 2000, 2005, 2007, 2008, 2010]) #LN yrs

#helper function to deal with nans
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/nan_helper.py')

if cropNam == 'maize2':
    months = [13,14,15,16] # mar, apr, may #[10,11,12,13,14,15]#
    states = ['PR','SP','MT','MS','GO','RS','SC','MG'] #'RS','SC','MG','RJ','ES'   
    saveNotes=''
    cropNam='maize';notes='2'
elif cropNam == 'maize1':
    months = [10,11,12,13,14,15] #jan, feb, mar, apr #[10,11,12,13,14,15]#
    states = ['PR','SP','RS','SC','GO','MT','MS','MG'] #'MT','MS','MG','RJ','ES'  
    saveNotes=''
    cropNam='maize';notes='1'
elif cropNam == 'maize':
    months = [10,11,12,13,14,15]#[12,13,14] #jan, feb, mar, apr
    states = ['PR','SP','RS','SC','GO','MT','MS','MG'] #'MT','MS','MG','RJ','ES'  
    notes='';saveNotes=''
elif cropNam == 'wheat':
    months = [5,6,7,8,9,10] #sep, Oct, Nov
    states = ['PR','SP','RS','SC','GO','MT','MS','MG']
    notes='';saveNotes=''
elif cropNam == 'soySouth':
    months = [-3,-2,-1,0,1,2]#[0,1] #ONDJFM
    states = ['PR','SP','RS','SC','GO','MT','MS','MG']#['PR','RS','SC']#,'GO','MS','MT']
    notes='';saveNotes='South';cropNam='soy'
    rustYrs = []#'2002','2003','2006']
elif cropNam == 'soyCerr':
    months = [-4,-3,-2,-1,0]#[-2,-1] #previous Nov, Dec
    states = ['PR','SP','RS','SC','GO','MT','MS','MG']
    notes='';saveNotes='Cerr';cropNam='soy'
    rustYrs = []#'2002','2003','2006']

if SMcorr =='T':
    SM_filter = 'SMfilter'
if SMcorr =='F':
    SM_filter = ''
    
monName= ['Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec']
 
s1 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/south1'+lag+'.npy')
s2 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/south2'+lag+'.npy')
s3 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/south3'+lag+'.npy')
cerr = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/cerrado'+lag+'.npy')
s0 = np.load('/Volumes/Data_Archive/Results/soilWaterBalance/reference'+lag+'.npy')
if smData is 'Dunne':
    fpath = '/Volumes/Data_Archive/Data/SoilMoisture/Dunne_soil/dunne_soil_1deg.txt'
    swc =  np.genfromtxt(fpath, dtype='f4') 
    swc = swc*10. #convert to mm
elif smData is 'IGBP':
    fpath = '/Volumes/Data_Archive/Data/SoilMoisture/IGBP-DIS/pawc_1deg.txt'
    swc =  np.genfromtxt(fpath, dtype='f4') 
smNC = netCDF4.Dataset('/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/raw/SoilM.nc')
sm = smNC.variables['SoilM'][smNC.variables['year'][:]>=1950,:,0:3,...] #start at 1950
sm = np.swapaxes(sm,0,1);sm=np.reshape(sm,[sm.shape[0]*sm.shape[1],sm.shape[2],sm.shape[3],sm.shape[4]],order='F')


sm2_DF = np.copy(sm[:,1,::-1,:])
smRZ_DF = (sm[:,0,::-1,:]+sm[:,1,::-1,:]+sm[:,2,::-1,:])
smMid_DF = (sm[:,0,::-1,:]+sm[:,1,::-1,:])

#read in major / minor cropping areas
if area == 'major':
    CCC = np.load('/Volumes/Data_Archive/Results/AmENSO/simpleCropMaps/halfPercMask/gpcc_grid/'+cropNam+'Map1.npy')
elif area == 'allCrop':
    CCC = np.load('/Volumes/Data_Archive/Results/AmENSO/simpleCropMaps/AllArea/gpcc_grid/'+cropNam+'Map1.npy')
elif area == 'all':
    CCC = np.load('/Volumes/Data_Archive/Results/AmENSO/simpleCropMaps/halfPercMask/gpcc_grid/'+cropNam+'Map1.npy')    
    CCC = np.ones([CCC.shape[0],CCC.shape[1]])
#elif area == 'MODIScrop':
    #placeholder for modis
CCC[CCC==13]=np.nan #set minor areas (code 13) equal to 0
CCC[CCC>0]=1
    
index = np.array(range(s1.shape[0]/12))*12
SMmonths = [10,11,12,13,14,15] #NDJFMA

fpath = '/Volumes/Data_Archive/Data/Precip/GPCC_nc/raw/'
f = netCDF4.Dataset(str(fpath+'sumP.nc'),'r',format='NETCDF4')

all_ind = np.array(range(years[0]-f.variables['year'][0],years[1]-f.variables['year'][0]))
en_ind = np.array(ENyrs-f.variables['year'][0]);ln_ind = np.array(LNyrs-f.variables['year'][0])
lats = f.variables['latitude'][:];lons = f.variables['longitude'][:]
lat_ind = np.where((lats<latMin)&(lats>latMax))[0]
lon_ind = np.where((lons<lonMax)&(lons>lonMin))[0]

p_all_0 = f.variables['sumP'][all_ind,:,lat_ind,lon_ind]
p_all_0 = np.nanmean(np.nanmean(np.nanmean(p_all_0,0),1),1)
p_all_1 = f.variables['sumP'][all_ind+1,:,lat_ind,lon_ind] #b/c the growing season straddles a year
p_all_1 = np.nanmean(np.nanmean(np.nanmean(p_all_1,0),1),1)
p_en_0 = f.variables['sumP'][en_ind,:,lat_ind,lon_ind]
p_en_0 = np.nanmean(np.nanmean(np.nanmean(p_en_0,0),1),1)
p_en_1 = f.variables['sumP'][en_ind+1,:,lat_ind,lon_ind] #b/c the growing season straddles a year
p_en_1 = np.nanmean(np.nanmean(np.nanmean(p_en_1,0),1),1)
p_ln_0 = f.variables['sumP'][ln_ind,:,lat_ind,lon_ind]
p_ln_0 = np.nanmean(np.nanmean(np.nanmean(p_ln_0,0),1),1)
p_ln_1 = f.variables['sumP'][ln_ind+1,:,lat_ind,lon_ind] #b/c the growing season straddles a year
p_ln_1 = np.nanmean(np.nanmean(np.nanmean(p_ln_1,0),1),1)
p_all = np.append(p_all_0,p_all_1,0)
p_en = np.append(p_en_0,p_en_1,0)
p_ln = np.append(p_ln_0,p_ln_1,0)

lon0 = 0;
lats = f.variables['latitude'][:] +.5; lons = f.variables['longitude'][:] -.5



if SMcorr == 'T':
    #get lats/lons
    fpath = '/Volumes/Data_Archive/Data/Precip/GPCC_nc/raw/'
    f = netCDF4.Dataset(str(fpath+'sumP.nc'),'r',format='NETCDF4')
    lats = f.variables['latitude'][:] + .5
    lons = f.variables['longitude'][:] - .5
    lons,lats = np.meshgrid(lons,lats) 
    ur_lat2 = 15.;ll_lat2 = -60;ll_lon2 = -85.;ur_lon2 = -30.
        #Define the monthly cropping coefficients (Kc) for sequences of 3 years
    Kc_s1 = np.array([np.nan,1.15,1.15,.5,0.1,np.nan,1.,1.,.5,0.1,np.nan,np.nan,
                 1.15,1.15,.5,0.1,np.nan,np.nan,1.,1.,0.1,0.1,np.nan,np.nan,
                 1.2,1.2,np.nan,.5,0.1,np.nan,np.nan,np.nan,1.15,1.15,0.35,0.1])
    Kc_s2 =   np.array([1.15,1.15,0.5,0.1,np.nan,np.nan,1.,1.,0.1,0.1,np.nan,np.nan,
                  1.2,1.2,np.nan,.5,0.1,np.nan,np.nan,np.nan,1.15,1.15,.35,0.1,
                  np.nan,1.15,1.15,0.5,0.1,np.nan,1.,1.,.5,0.1,np.nan,np.nan])
    Kc_s3 =  np.array([1.2,1.2,np.nan,0.5,0.1,np.nan,np.nan,np.nan,1.15,1.15,.35,0.1,
                 np.nan,1.15,1.15,0.5,0.1,np.nan,1.,1.,.5,0.1,np.nan,np.nan,
                 1.15,1.15,0.5,0.1,np.nan,np.nan,1.,1.,0.1,0.1,np.nan,np.nan])
    Kc_cer =  np.array([0.5,0.1,1.2,1.2,0.5,0.1,1.,0.1,0.5,np.nan,1.15,1.15,
                        0.5,0.1,1.2,1.2,0.5,0.1,1.,0.1,0.5,np.nan,1.15,1.15,
                        0.5,0.1,1.2,1.2,0.5,0.1,1.,0.1,0.5,np.nan,1.15,1.15])
    
    #tile to the full size of the dataset
    Kc_s1=np.tile(Kc_s1,3);Kc_s2=np.tile(Kc_s2,3)
    Kc_s3=np.tile(Kc_s3,3);Kc_cer=np.tile(Kc_cer,3)
    
    #linearly interpolate between months in the Kc values
    nans, x= nan_helper(Kc_s1)
    Kc_s1[nans] = np.interp(x(nans),x(~nans),Kc_s1[~nans])
    nans, x= nan_helper(Kc_s2)
    Kc_s2[nans] = np.interp(x(nans),x(~nans),Kc_s2[~nans])   
    nans, x= nan_helper(Kc_s3)
    Kc_s3[nans] = np.interp(x(nans),x(~nans),Kc_s3[~nans])
    nans, x= nan_helper(Kc_cer)
    Kc_cer[nans] = np.interp(x(nans),x(~nans),Kc_cer[~nans])
        
    Kc_s = (Kc_s1+Kc_s2+Kc_s3)/3
    
    smMidAnom = np.zeros([index.size*np.size(SMmonths),smMid_DF.shape[1],smMid_DF.shape[2]])*np.nan
    sm2Anom = np.zeros([index.size*np.size(SMmonths),sm2_DF.shape[1],sm2_DF.shape[2]])*np.nan
    smRZanom = np.zeros([index.size*np.size(SMmonths),smRZ_DF.shape[1],smRZ_DF.shape[2]])*np.nan
    s0anom = np.zeros([index.size*np.size(SMmonths),s0.shape[1],s0.shape[2]])*np.nan
    for i in range(np.size(SMmonths)):
        im = SMmonths[i]-1
        if(im<0):
            indMon = index+im+12
        elif(im>11):
            indMon = index+im-12
        else:
            indMon = index+im
        ind = np.arange(i,s0anom.shape[0],np.size(SMmonths))
        smRZanom[ind,...] = smRZ_DF[indMon,...] - np.nanmean(smRZ_DF[indMon,...],axis=0)
        smMidAnom[ind,...] = smMid_DF[indMon,...] - np.nanmean(smMid_DF[indMon,...],axis=0)
        sm2Anom[ind,...] = sm2_DF[indMon,...] - np.nanmean(sm2_DF[indMon,...],axis=0)
        s0anom[ind,...] = s0[indMon,...] - np.nanmean(s0[indMon,...],axis=0)
    
    smRZcorr = np.zeros(np.shape(swc))*np.nan
    smMidCorr = np.zeros(np.shape(swc))*np.nan
    sm2Corr = np.zeros(np.shape(swc))*np.nan
    for row in range(smRZ_DF.shape[1]):
        for col in range(smRZ_DF.shape[2]):
            smRZcorr[row,col]=np.corrcoef(smRZanom[:,row,col],s0anom[:,row,col])[1][0]
            smMidCorr[row,col]=np.corrcoef(smMidAnom[:,row,col],s0anom[:,row,col])[1][0]
            sm2Corr[row,col]=np.corrcoef(sm2Anom[:,row,col],s0anom[:,row,col])[1][0]
    smMidCorr = np.ma.masked_invalid(smMidCorr)  


    monNums = np.arange(0,14,2)
    kcNums = np.arange(0.2,1.2,.4)
    monNams= ['Sep','Nov','Jan','Mar','May','July','Sep']
    
    fig = plt.figure();
    ax1=plt.subplot(122);ax2=plt.subplot(521);ax3=plt.subplot(523);ax4=plt.subplot(525);
    ax5=plt.subplot(529);ax6=plt.subplot(527)
    m1 = Basemap(projection='mill',lon_0=-0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1); 
    m1.drawcoastlines();m1.drawcountries();
    im1 = m1.pcolormesh(lons,lats,smMidCorr,latlon=True,shading='flat',vmin=-1,vmax=1,cmap='PiYG');
    ax6.bar(range(p_en.size),p_en-p_all,width=0.33,color='darkred');
    ax6.bar(np.array(range(p_ln.size))+0.33,p_ln-p_all,width=0.33,color='cornflowerblue');
    ax6.set_xlim(8,20);ax6.set_ylim(-20,35)
    ax2.plot(Kc_s1[8:21],'-k',lw=2);ax3.plot(Kc_s2[8:21],'-k',lw=2);ax4.plot(Kc_s3[8:21],'-k',lw=2)
    ax5.plot(Kc_cer[8:21],'k',lw=2) ;ax5.set_yticks(kcNums);ax2.set_xticklabels('')
    ax2.set_yticks(kcNums);ax3.set_yticks(kcNums);ax4.set_yticks(kcNums)
    ax5.set_xticks(monNums);ax5.set_xticklabels(monNams);
    ax4.set_xticks(monNums);ax4.set_xticklabels(monNams);
    cb1 = m1.colorbar(im1,"right", size="5%", pad="2%");
    ax3.set_ylim([0,1.3]);ax2.set_ylim([0,1.3]);ax4.set_ylim([0,1.3]);ax5.set_ylim([0,1.3])
    ax3.set_xlim([0,12]);ax2.set_xlim([0,12]);ax4.set_xticklabels('');ax6.set_xticklabels('')
    ax5.set_xlim([0,12]);ax4.set_xlim([0,12]);ax3.set_xticklabels('')
    ax3.set_ylabel('Traditional cropping cycle Kc curve\n\n year 2')
    ax2.set_ylabel('year 1');ax4.set_ylabel('year 3')
    ax5.set_ylabel('Safrinha cycle\n Kc curve\n')
    ax6.set_ylabel('Precip.\nanomalies\n')
    ax1.set_title('Soil water content - soil moisture \n correlation coefficient')
    fig.set_size_inches(12,8)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/Kc_curves_pAnom.png')
    plt.close()
    
    fig = plt.figure();
    ax1=plt.subplot(122);ax2=plt.subplot(521);ax3=plt.subplot(523);ax4=plt.subplot(525);
    ax5=plt.subplot(529);ax6=plt.subplot(527)
    m1 = Basemap(projection='mill',lon_0=-0,lat_0=0,llcrnrlat=ll_lat2,urcrnrlat=ur_lat2,llcrnrlon=ll_lon2,urcrnrlon=ur_lon2,ax=ax1); 
    m1.drawcoastlines();m1.drawcountries();
    im1 = m1.pcolormesh(lons,lats,smMidCorr,latlon=True,shading='flat',vmin=-1,vmax=1,cmap='PiYG');
    ax6.plot(p_all,'--k',lw=2);
    ax6.plot(p_ln,lw=3,color='cornflowerblue');ax6.plot(p_en,lw=3,color='darkred');
    ax6.set_xlim(8,20);ax6.set_ylim(25,260)
    ax2.plot(Kc_s1[8:21],'-k',lw=2);ax3.plot(Kc_s2[8:21],'-k',lw=2);ax4.plot(Kc_s3[8:21],'-k',lw=2)
    ax5.plot(Kc_cer[8:21],'k',lw=2) ;ax5.set_yticks(kcNums);ax2.set_xticklabels('')
    ax2.set_yticks(kcNums);ax3.set_yticks(kcNums);ax4.set_yticks(kcNums)
    ax5.set_xticks(monNums);ax5.set_xticklabels(monNams);
    ax4.set_xticks(monNums);ax4.set_xticklabels(monNams);
    cb1 = m1.colorbar(im1,"right", size="5%", pad="2%");
    ax3.set_ylim([0,1.3]);ax2.set_ylim([0,1.3]);ax4.set_ylim([0,1.3]);ax5.set_ylim([0,1.3])
    ax3.set_xlim([0,12]);ax2.set_xlim([0,12]);ax4.set_xticklabels('');ax6.set_xticklabels('')
    ax5.set_xlim([0,12]);ax4.set_xlim([0,12]);ax3.set_xticklabels('')
    ax3.set_ylabel('Traditional cropping cycle Kc curve\n\n year 2')
    ax2.set_ylabel('year 1');ax4.set_ylabel('year 3')
    ax5.set_ylabel('Safrinha cycle\n Kc curve\n')
    ax6.set_ylabel('Precipitation\n')
    ax1.set_title('Soil water content - soil moisture \n correlation coefficient')
    fig.set_size_inches(12,8)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/SWC/Kc_curves_pAbs.png')
    plt.close()
else:
    smRZcorr = np.ones(np.shape(swc))
    smMidCorr = np.ones(np.shape(swc)) 
    sm2Corr = np.ones(np.shape(swc)) 


if SMperc == 'T':
    smMid_DF = smMid_DF / swc[:-30,:]
    smRZ_DF = smMid_DF / swc[:-30,:]
    sm2_DF = sm2_DF / swc[:-30,:]

#Don't include the filter by soil moisture, this seems too much messing
smMidCorr[smMidCorr>=0.5]=1; smMidCorr[smMidCorr<0.5]=0
smRZcorr[smRZcorr>=0.5]=1; smRZcorr[smRZcorr<0.5]=0

#Convertt to percent anomalies
s1 = s1*CCC*smMidCorr/swc; s2 = s2*CCC*smMidCorr/swc; s3 = s3*CCC*smMidCorr/swc; cerr = cerr*CCC*smMidCorr/swc


crop = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/AmENSO_yld/AmENSO_BR_'+cropNam+notes+'_allGDD2allstate.csv')
crop = crop.set_index('state',drop=True,append=True)

sheffMircaCodes = np.load('/Volumes/Data_Archive/Data/CropCalendar/'+
                    'CCCmask/Sheffield_grid/UnitCode/SubnatUnitCode.npy')
adminCodes = pd.read_table('/Volumes/Data_Archive/Data/CropCalendar/unit_name_state.txt',sep='\t')
#re-parse the unit code names
for n in range(0,np.size(adminCodes['UnitName'])):
    if np.size(adminCodes['UnitName'][n].split('_'))>1:
	newName = adminCodes['UnitName'][n].split('_')[1]
	adminCodes['UnitName'][n]=newName

mismatch = []

s0Mean = []
sMean_st =  []
sMean_SWC = []
sMean_Yld = []
cer_SWC = []

if (any(mo>11 for mo in months)|any(mo<0 for mo in months)):
    s0AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    smMidAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    smRZAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    s1AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    s2AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    s3AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    sMeanAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    cerrAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
else:
    s0AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
    smMidAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
    smRZAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
    s1AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
    s2AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
    s3AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
    sMeanAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
    cerrAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan

for state in states:
    sheffMircaCodesTemp = copy.deepcopy(sheffMircaCodes)
    code = adminCodes['UnitCode'][adminCodes['UnitName']==state.title()]
    if code.size==0:
        mismatch.append(state)
        continue
    for num in code:
        sheffMircaCodesTemp[sheffMircaCodes==num]=1
        sheffMircaCodesSM = sheffMircaCodesTemp[:-30,:]
    sheffMircaCodesTemp[sheffMircaCodesTemp!=1]=0   
    sheffMircaCodesSM[sheffMircaCodesSM!=1]=0  
    print mismatch
    for im in range(0,np.size(months)):
        yields = crop[yldAnom].loc[:,state]
        newInd = index+months[im]
        if any(mo>11 for mo in months):
            yields = yields['1950':'2009']
            newInd=newInd[:-1]
        elif any(mo<0 for mo in months):
            yields = yields['1951':'2010']
            newInd=newInd[1:]
        else:
            yields = yields['1950':'2010']
        s0st = np.nanmean(s0[:,sheffMircaCodesTemp==1],axis=1)
        smMid = np.nanmean(smMid_DF[:,sheffMircaCodesSM==1],axis=1)
        smRZ = np.nanmean(smRZ_DF[:,sheffMircaCodesSM==1],axis=1)
        ws1 = np.nanmean(s1[:,sheffMircaCodesTemp==1],axis=1)
        ws2 =  np.nanmean(s2[:,sheffMircaCodesTemp==1],axis=1)
        ws3 =  np.nanmean(s3[:,sheffMircaCodesTemp==1],axis=1)
        wsMean = np.nanmean([ws1,ws2,ws3],0)
        wcerr =  np.nanmean(cerr[:,sheffMircaCodesTemp==1],axis=1)
        s0st=s0st[newInd]
        ws1 = ws1[newInd];ws2 = ws2[newInd];ws3 = ws3[newInd]
        wcerr = wcerr[newInd];wsMean = wsMean[newInd];       
        smMid=smMid[newInd];smRZ=smRZ[newInd]
        s0AllMon[im,...]=s0st
        smMidAllMon[im,...]=smMid
        smRZAllMon[im,...]=smRZ
        s1AllMon[im,...] = ws1
        s2AllMon[im,...] = ws2
        s3AllMon[im,...] = ws3
        sMeanAllMon[im,...] = wsMean
        cerrAllMon[im,...] = wcerr
        #if the crop is soy, take out the rust years
        if cropNam == 'soy':
            for yr in rustYrs:
                yields[yr] = np.nan
        #limit to non-nan
        nonNan = np.array(~np.isnan(yields))
        smMid=smMid[nonNan]
        smRZ=smRZ[nonNan]
        yields=yields[nonNan]        
        ws1=ws1[nonNan]
        ws2=ws2[nonNan]
        ws3=ws3[nonNan]
        wsMean=wsMean[nonNan]
        wcerr=wcerr[nonNan]
        print monName[months[im]%12]
        print 's1 '+state+'   '+np.str(np.corrcoef(ws1,yields)[0][1])
        print 's2 '+state+'   '+np.str(np.corrcoef(ws2,yields)[0][1])
        print 's3 '+state+'   '+np.str(np.corrcoef(ws3,yields)[0][1])
        print 'cerr '+state+'   '+np.str(np.corrcoef(wcerr,yields)[0][1])
    s0AllMon=s0AllMon[:,nonNan]
    smMidAllMon=smMidAllMon[:,nonNan]
    smRZAllMon=smRZAllMon[:,nonNan]
    s1AllMon=s1AllMon[:,nonNan]
    s2AllMon=s2AllMon[:,nonNan]
    s3AllMon=s3AllMon[:,nonNan]
    sMeanAllMon=sMeanAllMon[:,nonNan]
    cerrAllMon=cerrAllMon[:,nonNan]
    s0_seasMean = np.nanmean(s0AllMon,axis=0)
    smMid_seasMean = np.nanmean(smMidAllMon,axis=0)
    smRZ_seasMean = np.nanmean(smRZAllMon,axis=0)
    s1_seasMean = np.nanmean(s1AllMon,axis=0)
    s2_seasMean = np.nanmean(s2AllMon,axis=0)
    s3_seasMean = np.nanmean(s3AllMon,axis=0) 
    sMean_seasMean = np.nanmean(sMeanAllMon,axis=0) 
    cerr_seasMean = np.nanmean(cerrAllMon,axis=0)

    sMean_SWC.append(sMean_seasMean)    
    cer_SWC.append(cerr_seasMean)
    sMean_Yld.append(yields.values)
    sMean_st.append(np.repeat(state,sMean_seasMean.shape[0]))  
    print 'Season Mean'
#    print 's1 '+state+'   '+np.str(np.corrcoef(s1_seasMean,yields)[0][1])
#    print 's2 '+state+'   '+np.str(np.corrcoef(s2_seasMean,yields)[0][1])
#    print 's3 '+state+'   '+np.str(np.corrcoef(s3_seasMean,yields)[0][1])
    print 'sMean '+state+'   '+np.str(np.corrcoef(sMean_seasMean,yields)[0][1])
    print 'cerr '+state+'   '+np.str(np.corrcoef(cerr_seasMean,yields)[0][1])
    print 'ET0 '+state+'   '+np.str(np.corrcoef(s0_seasMean,smMid_seasMean)[0][1])
    print 'soilM '+state+'   '+np.str(np.corrcoef(yields,smMid_seasMean)[0][1])
    print ' '
    fig = plt.figure();
    plt.plot(yields.index.year,s1_seasMean-np.nanmean(s1_seasMean),'--k')
    plt.plot(yields.index.year,s2_seasMean-np.nanmean(s2_seasMean),'--k')
    plt.plot(yields.index.year,s3_seasMean-np.nanmean(s3_seasMean),'--k')
    plt.plot(yields.index.year,sMean_seasMean-np.nanmean(sMean_seasMean),'--k',lw=2)
    plt.plot(yields.index.year,cerr_seasMean-np.nanmean(cerr_seasMean),'--k',lw=2)
    plt.plot(yields.index.year,yields.values/10,'-o',color='orange',lw=2)
    plt.plot(yields.index.year,(smRZ_seasMean-np.nanmean(smRZ_seasMean)),'b',lw=2)
    #(sMean_seasMean,cerr_seasMean,yields.values),('south','cerrado','yields'))
    fig.set_size_inches(10, 6)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/safrinhaMaize/SWC_TS/'+cropNam+notes+saveNotes+state+yldAnom+'_swc'+lag+area+'.png')
    plt.close()
    
    fig = plt.figure();
    plt.plot(yields.index.year,s1_seasMean-np.nanmean(s1_seasMean),'--k')
    plt.plot(yields.index.year,s2_seasMean-np.nanmean(s2_seasMean),'--k')
    plt.plot(yields.index.year,s3_seasMean-np.nanmean(s3_seasMean),'--k')
    plt.plot(yields.index.year,sMean_seasMean-np.nanmean(sMean_seasMean),'--k',lw=2)
    plt.plot(yields.index.year,cerr_seasMean-np.nanmean(cerr_seasMean),'darkgreen',lw=2)
    plt.plot(yields.index.year,(smRZ_seasMean-np.nanmean(smRZ_seasMean)),'b',lw=2)
    fig.set_size_inches(10, 6)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/safrinhaMaize/SWC_TS/'+cropNam+notes+saveNotes+state+'sm_swc'+lag+area+'.png')
    plt.close()

    fig = plt.figure();
#    plt.plot(yields.index.year,s1_seasMean-np.nanmean(s1_seasMean),'--k')
#    plt.plot(yields.index.year,s2_seasMean-np.nanmean(s2_seasMean),'--k')
#    plt.plot(yields.index.year,s3_seasMean-np.nanmean(s3_seasMean),'--k')
#    plt.plot(yields.index.year,sMean_seasMean-np.nanmean(sMean_seasMean),'--k',lw=2)
    plt.plot(yields.index.year,(cerr_seasMean-np.nanmean(cerr_seasMean))/np.std(cerr_seasMean),'--k',lw=2)
    plt.plot(yields.index.year,(yields.values-np.mean(yields.values))/np.std(yields.values),'-o',color='orange',lw=2)
    plt.plot(yields.index.year,(smMid_seasMean-np.nanmean(smMid_seasMean))/np.std(smMid_seasMean),'b',lw=2)
    #(sMean_seasMean,cerr_seasMean,yields.values),('south','cerrado','yields'))
    fig.set_size_inches(10, 6)
    fig.savefig('/Users/westonanderson/Desktop/Columbia/Research/Results/AmENSO_yields/safrinhaMaize/SWC_TS/'+cropNam+notes+saveNotes+state+'_SM_Yld_SWC'+lag+area+'.png')
    plt.close()
    
    if (any(mo>11 for mo in months)|any(mo<0 for mo in months)):
        s0AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
        smMidAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
        smRZAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
        s1AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
        s2AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
        s3AllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
        sMeanAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
        cerrAllMon = np.zeros([np.size(months),index.shape[0]-1])*np.nan
    else:
        s0AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
        smMidAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
        smRZAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
        s1AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
        s2AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
        s3AllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
        sMeanAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
        cerrAllMon = np.zeros([np.size(months),index.shape[0]])*np.nan
 
#plt.scatter(np.ravel(sMean_SWC),np.ravel(sMean_Yld))
#plt.text(0.7,0.2,np.corrcoef(np.ravel(sMean_SWC),np.ravel(sMean_Yld))[0][1]) 
#plt.scatter(np.ravel(cer_SWC),np.ravel(sMean_Yld))
#plt.text(0.7,0.2,np.corrcoef(np.ravel(cer_SWC),np.ravel(sMean_Yld))[0][1])   
#sMean_allSt = pd.DataFrame(columns = ['state','SWC','YieldAnom']) 
