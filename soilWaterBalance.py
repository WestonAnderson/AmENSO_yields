# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:19:06 2016

Code to calculate the soil water balance. Kc curves defined for Brazilian cropping systems

Input data:
GLDAS soil moisture (Noah) to initialize the balance
Sheffield climate data
Dunne et al. (or IGBP) Plant available water capacity in the soil
PET was calculated from the sheffield data and GLDAS data following FAO56

To generalize this code, make if for one Kc and take out the repetition

@author: westonanderson

EDIT: 09/12/16 - removed the plotting portions of the script and made them two new scripts

EDIT: 11/16/16 - Added a check to make sure actual soil moisture (in addition to soil moisture at last step)
                never exceeds SWC. This didn't change results for anything (b/c last SM object linked to current one in python)
"""
import netCDF4
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
#plt.ioff()

anom = '' #'anom' or ''
percent = 'pct' #'pct' or ''
lag = 0 #lag or ''
smData = 'IGBP' #Dunne or IGBP
years = [1950,2011]
numYrs = years[1]-years[0]


#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#
#read in a function to regrid from one grid to another
# function is: ### regrid_ASCII(oldGrid,degOld,lon0,latMin,latMax,degNew,sumORavg)
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/nan_helper.py')
execfile('/Users/westonanderson/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/animate_array.py')
#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#

#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#
#read in the plant available soil water content
if smData is 'Dunne':
    fpath = '/Volumes/Data_Archive/Data/SoilMoisture/Dunne_soil/dunne_soil_1deg.txt'
    swc =  np.genfromtxt(fpath, dtype='f4') 
    swc = swc*10. #convert to mm
elif smData is 'IGBP':
    fpath = '/Volumes/Data_Archive/Data/SoilMoisture/IGBP-DIS/pawc_1deg.txt'
    swc =  np.genfromtxt(fpath, dtype='f4') 

#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#
# Read in the PET data (ET0)
fpath =  '/Volumes/Data_Archive/Results/AmENSO/IntegratedVars/monthlyGlobal/'
ncdfVar = netCDF4.Dataset(str(fpath+'PET.nc'),'r',format='NETCDF4')
ET0 = ncdfVar.variables['PET'][years[0]-ncdfVar.variables['year'][0]:,:,...]#starts in 1948, start in 1950 
#reorder ET to be a single time axis of continuous months
ET0 = np.swapaxes(ET0,0,1);ET0 = np.reshape(ET0,[ET0.shape[0]*ET0.shape[1],ET0.shape[2],ET0.shape[3]],order='F')
#adjust the ET longitudes to match the SWC
ET0 = np.append(ET0[:,::-1,180:],ET0[:,::-1,:180],axis=2)
#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#
## Read in the precipitation data
#fpath = '/Volumes/Data_Archive/Data/Sheffield/'
#f = netCDF4.Dataset(str(fpath+'prcp_monthly_1948-2010.nc'),'r',format='NETCDF4')
#p = f.variables['prcp'][24:,...]#starts in 1948, start in 1950 
##convert from kg/(s*m^2) to mm/month. But it looks like the units are wrong. Need to convert .1 kg/s to mm/month
#days = [31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.,
# 31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.,
# 31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.,
# 31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.,] #days in months. 1952 was a leap year
#days = np.tile(days,16);days=days[:np.shape(p)[0]]
#p = np.squeeze(p)
#p =  np.append(p[:,::-1,180:],p[:,::-1,:180],axis=2)
#for ir in range(p.shape[1]):
#    for ic in range(p.shape[2]):
#         p[:,ir,ic]=p[:,ir,ic]*days*86400.*10
##save lat and lon for later
#lats = f.variables['latitude'][:]
#lats = lats[::-1]
#lons = f.variables['longitude'][:]
fpath = '/Volumes/Data_Archive/Data/Precip/GPCC_nc/raw/'
f = netCDF4.Dataset(str(fpath+'sumP.nc'),'r',format='NETCDF4')
p = f.variables['sumP'][(years[0]-f.variables['year'][0]):(numYrs+(years[0]-f.variables['year'][0])),...]#start at 1950
p = np.swapaxes(p,0,1);p= np.reshape(p,[p.shape[0]*p.shape[1],p.shape[2],p.shape[3]],order='F')
p[p.mask]=0
lats = f.variables['latitude'][:] +.5
lons = f.variables['longitude'][:] -.5
lon0 = 0;
#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#
# Read in the soil moisture data to use to initialize the soil water balance
fpath = '/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/raw/'
f = netCDF4.Dataset(str(fpath+'SoilM.nc'),'r',format='NETCDF4')
sm = f.variables['SoilM'][0,years[0]-f.variables['year'][0],0:3,...] #Jan soil moisture from 0-1m
sm = np.nansum(sm,axis=0)*0.5 #estimate plant available soil moisture as 1/2 of column integrated soil moisture
sm = sm[::-1,:]
sm0 = np.zeros([180,360])*np.nan
sm0[:-30,:]=sm


months =range(36)
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
Kc_s1=np.tile(Kc_s1,(numYrs/3)+3);Kc_s2=np.tile(Kc_s2,(numYrs/3)+3)
Kc_s3=np.tile(Kc_s3,(numYrs/3)+3);Kc_cer=np.tile(Kc_cer,(numYrs/3)+3)

#linearly interpolate between months in the Kc values
nans, x= nan_helper(Kc_s1)
Kc_s1[nans] = np.interp(x(nans),x(~nans),Kc_s1[~nans])
nans, x= nan_helper(Kc_s2)
Kc_s2[nans] = np.interp(x(nans),x(~nans),Kc_s2[~nans])   
nans, x= nan_helper(Kc_s3)
Kc_s3[nans] = np.interp(x(nans),x(~nans),Kc_s3[~nans])
nans, x= nan_helper(Kc_cer)
Kc_cer[nans] = np.interp(x(nans),x(~nans),Kc_cer[~nans])

Kc_s1 = Kc_s1[36:(3+numYrs)*12]; Kc_s2 = Kc_s2[36:(3+numYrs)*12]
Kc_s3 = Kc_s3[36:(3+numYrs)*12]; Kc_cer = Kc_cer[36:(3+numYrs)*12]


Kc_s1 = np.roll(Kc_s1,lag)
Kc_s2 = np.roll(Kc_s2,lag)
Kc_s3 = np.roll(Kc_s3,lag)
Kc_cer = np.roll(Kc_cer,lag)

#make everything 1d so that boolean indexing is straightforward
sm0 = np.reshape(sm0,[sm0.shape[0]*sm0.shape[1]])
swc = np.reshape(swc,[swc.shape[0]*swc.shape[1]])
p = np.reshape(p,[p.shape[0],p.shape[1]*p.shape[2]])
ET0 = np.reshape(ET0,[ET0.shape[0],ET0.shape[1]*ET0.shape[2]])
w0 = np.zeros([12*numYrs,180*360])*np.nan
w_s1 = np.zeros([12*numYrs,180*360])*np.nan
w_s2 = np.zeros([12*numYrs,180*360])*np.nan
w_s3 = np.zeros([12*numYrs,180*360])*np.nan
w_cer = np.zeros([12*numYrs,180*360])*np.nan
E0 = np.zeros([12*numYrs,180*360])*np.nan
Es1 = np.zeros([12*numYrs,180*360])*np.nan
Es2 = np.zeros([12*numYrs,180*360])*np.nan
Es3 = np.zeros([12*numYrs,180*360])*np.nan
Ecer = np.zeros([12*numYrs,180*360])*np.nan
Rcer = np.zeros([12*numYrs,180*360])*np.nan
Rs1 =np.zeros([12*numYrs,180*360])*np.nan
Rs2 =np.zeros([12*numYrs,180*360])*np.nan
Rs3 =np.zeros([12*numYrs,180*360])*np.nan
swc[swc==0]=np.nan#in denom, so can't divide by 0
sm0[swc==np.nan]=np.nan
sm0[sm0>swc]=swc[sm0>swc] #soil can never hold more than swc


#start the soil water balance calculation
for i in range(0,numYrs*12):
    if i==0: #initialize the first step with the soil moisture estimate
        w_lst_0 =sm0        
        w_lst_s1=sm0
        w_lst_s2=sm0
        w_lst_s3=sm0
        w_lst_cer=sm0
    else: #define the last soil moisture state and runoff
        w_lst_0 = w0[i-1,...]
        w_lst_s1 = w_s1[i-1,...]
        w_lst_s2 = w_s2[i-1,...]
        w_lst_s3 = w_s3[i-1,...]
        w_lst_cer = w_cer[i-1,...]
        Rcer[i-1,...] = np.maximum(0,w_lst_cer-swc)
        Rs1[i-1,...] = np.maximum(0,w_lst_s1-swc)
        Rs2[i-1,...] = np.maximum(0,w_lst_s2-swc)
        Rs3[i-1,...] = np.maximum(0,w_lst_s3-swc)
        w_lst_0[w_lst_0>swc]=swc[w_lst_0>swc] #set max available soil moisture to the holding capacity (so account for runoff)
        w_lst_s1[w_lst_s1>swc]=swc[w_lst_s1>swc]        
        w_lst_s2[w_lst_s2>swc]=swc[w_lst_s2>swc]
        w_lst_s3[w_lst_s3>swc]=swc[w_lst_s3>swc]
        w_lst_cer[w_lst_cer>swc]=swc[w_lst_cer>swc]
        w_lst_0[w_lst_0<0]=0 #set min available soil moisture to 0
        w_lst_s1[w_lst_s1<0]=0
        w_lst_s2[w_lst_s2<0]=0 
        w_lst_s3[w_lst_s3<0]=0 
        w_lst_cer[w_lst_cer<0]=0 
        w0[i-1,w0[i-1,:]>swc]=swc[w0[i-1,:]>swc]
        w_s1[i-1,w_s1[i-1,:]>swc]=swc[w_s1[i-1,:]>swc] #set max available soil moisture to the holding capacity (so account for runoff)
        w_s2[i-1,w_s2[i-1,:]>swc]=swc[w_s2[i-1,:]>swc]
        w_s3[i-1,w_s3[i-1,:]>swc]=swc[w_s3[i-1,:]>swc]
        w_cer[i-1,w_cer[i-1,:]>swc]=swc[w_cer[i-1,:]>swc]
        w0[i-1,w0[i-1,:]<0]=0 #set min available soil moisture to 0
        w_s1[i-1,w_s1[i-1,:]<0]=0
        w_s2[i-1,w_s2[i-1,:]<0]=0 
        w_s3[i-1,w_s3[i-1,:]<0]=0 
        w_cer[i-1,w_cer[i-1,:]<0]=0 
    f0 = w_lst_0/swc #factor to determine the sm stress, which modifies ET0
    fs1 = w_lst_s1/swc
    fs2 = w_lst_s2/swc
    fs3 = w_lst_s3/swc
    fcer = w_lst_cer/swc
    E0[i,...] = f0*ET0[i,...]  #determine final ET based on PET (ET0), crop growth stage and sm stress
    Es1[i,...] = fs1*ET0[i,...]*Kc_s1[i]
    Es2[i,...] = fs2*ET0[i,...]*Kc_s2[i]
    Es3[i,...] = fs3*ET0[i,...]*Kc_s3[i]
    Ecer[i,...] = fcer*ET0[i,...]*Kc_cer[i]
    w_now_0 = p[i,...]- E0[i,...] + w_lst_0  #precip - evap + soil moisture
    w_now_s1 = p[i,...]- Es1[i,...] + w_lst_s1
    w_now_s2 = p[i,...]- Es2[i,...] + w_lst_s2
    w_now_s3 = p[i,...]- Es3[i,...] + w_lst_s3
    w_now_cer = p[i,...]- Ecer[i,...] + w_lst_cer
    w0[i,...]=w_now_0 #write to the available soil moisture dataframes
    w_s1[i,...]=w_now_s1
    w_s2[i,...]=w_now_s2
    w_s3[i,...]=w_now_s3
    w_cer[i,...]=w_now_cer
w_lst_0 = w0[i,...]
w_lst_s1 = w_s1[i,...]
w_lst_s2 = w_s2[i,...]
w_lst_s3 = w_s3[i,...]
w_lst_cer = w_cer[i,...]
#record the last values of Runoff
Rcer[i,...] = np.maximum(0,w_lst_cer-swc)
Rs1[i,...] = np.maximum(0,w_lst_s1-swc)
Rs2[i,...] = np.maximum(0,w_lst_s2-swc)
Rs3[i,...] = np.maximum(0,w_lst_s3-swc)
w_lst_s1[w_lst_s1>swc]=swc[w_lst_s1>swc] #set max available soil moisture to the holding capacity (so account for runoff)
w_lst_0[w_lst_0>swc]=swc[w_lst_0>swc]
w_lst_s2[w_lst_s2>swc]=swc[w_lst_s2>swc]
w_lst_s3[w_lst_s3>swc]=swc[w_lst_s3>swc]
w_lst_cer[w_lst_cer>swc]=swc[w_lst_cer>swc]
w_lst_0[w_lst_0<0]=0 #set min available soil moisture to 0
w_lst_s1[w_lst_s1<0]=0
w_lst_s2[w_lst_s2<0]=0 
w_lst_s3[w_lst_s3<0]=0 
w_lst_cer[w_lst_cer<0]=0 

ETs = (Es1+Es2+Es3)/3
Rs = (Rs1+Rs2+Rs3)/3
w_s = (w_s1+w_s2+w_s3)/3

#reshape the variables
w0 = np.reshape(w0,[w0.shape[0],180,360])
w_s1 = np.reshape(w_s1,[w_s1.shape[0],180,360])
w_s2 = np.reshape(w_s2,[w_s2.shape[0],180,360])
w_s3 = np.reshape(w_s3,[w_s3.shape[0],180,360])
w_cer = np.reshape(w_cer,[w_cer.shape[0],180,360])
Ecer = np.reshape(Ecer,[Ecer.shape[0],180,360])
ETs = np.reshape(ETs,[ETs.shape[0],180,360])
Rcer = np.reshape(Rcer,[Rcer.shape[0],180,360])
Rs = np.reshape(Rs,[Rs.shape[0],180,360])
ET0 = np.reshape(ET0,[ET0.shape[0],180,360])
swc = np.reshape(swc,[180,360])
p = np.reshape(p,[p.shape[0],180,360])

#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#~*~#
######## Calculate anomalies wrt the reference ET ###############
if anom=='anom':
    w_s1 =w_s1-w0
    w_s2 =w_s2-w0
    w_s3 =w_s3-w0
    w_cer =w_cer-w0
if percent=='pct':
    w_s1 =w_s1/swc
    w_s2 =w_s2/swc
    w_s3 =w_s3/swc
    w_cer =w_cer/swc

#write the objects out
np.save('/Volumes/Data_Archive/Results/soilWaterBalance/reference_lag'+str(lag)+anom+percent+'.npy',w0)
np.save('/Volumes/Data_Archive/Results/soilWaterBalance/south1_lag'+str(lag)+anom+percent+'.npy',w_s1)
np.save('/Volumes/Data_Archive/Results/soilWaterBalance/south2_lag'+str(lag)+anom+percent+'.npy',w_s2)
np.save('/Volumes/Data_Archive/Results/soilWaterBalance/south3_lag'+str(lag)+anom+percent+'.npy',w_s3)
np.save('/Volumes/Data_Archive/Results/soilWaterBalance/cerrado_lag'+str(lag)+anom+percent+'.npy',w_cer)

