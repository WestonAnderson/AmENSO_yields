# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 23:33:01 2016
@author: westonanderson

Estimates of ENSO induced yield variability trends

Analysis Structure
1. Establish that trends in production variability exist
2. Establish that ENSO is the dominant driver of production variability
3. Remove ENSO-induced changes in production, look at remaining trends in variability
4. Estimate possible trends in ENSO-induced production variability via
	4.0. Phase randomization of production PCs directly
      4.1  Regress onto current yields, phase randomize ENSO
	4.2. Regression onto paleo estimates of ENSO
	4.3. Regression onto CMIP-5 estimates of ENSO?
5. Compare observed production variability over the latter period (1996-2010) 
    and early period (1980-1995) to possible production variability in 15 yr windows


EDIT - 12/17/16 - script created to compare to the multi-EOF function

"""
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from scipy import stats
import matplotlib.cm as cm
start = time.clock()

pvalThresh = 0.1
crops = ['maize']#,'maize','soy','wheat']#

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                 Helper Functions                  #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/project/AmENSO_yields/corrST_'+str(pvalThresh)+'_dict.py').read())
#GADM dict to link to shapefiles
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_dict.py').read())
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#        Calculate trends in variability            #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
corrST_dict = {'wheat':wheatSTs,'maize':maizeSTs,'soy':soySTs}

clrBr = np.arange(-1,1.1,.1)
norm = Normalize(vmin=-1, vmax=1, clip=True)
mapper=cm.ScalarMappable(norm=norm, cmap='BrBG')

for crop in crops:
    corSTs = corrST_dict[crop]
    plt.figure()
    m = Basemap(llcrnrlon=-90,llcrnrlat=-60,urcrnrlon=-20,urcrnrlat=10,
            projection='merc',lon_0=0)#,lat_1=33,lat_2=45,lon_0=-95)
    BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                    name='states', drawbounds=True)
    # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if shape_dict['HASC_1'] in corSTs.keys():
            colors[shape_dict['HASC_1']]=corSTs[shape_dict['HASC_1']]/2+0.5
        else:colors[shape_dict['HASC_1']]=0.5
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        poly = Polygon(seg,facecolor=cm.BrBG(colors[names[nshape]]), edgecolor='k')
        ax.add_patch(poly) 
    ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                    name='states', drawbounds=True)
    # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if shape_dict['HASC_1'] in corSTs.keys():
            colors[shape_dict['HASC_1']]=corSTs[shape_dict['HASC_1']]/2+0.5
        else:colors[shape_dict['HASC_1']]=0.5
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        poly = Polygon(seg,facecolor=cm.BrBG(colors[names[nshape]]), edgecolor='k')
        ax.add_patch(poly) 
    mapper.set_array(clrBr);plt.colorbar(mapper)
    plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/AmENSO_yields/Yields/yieldCorr/corrMaps/'+crop+'/SouthAmerica_'+str(pvalThresh)+'.png')
    plt.close() 
    
    
    if crop is 'maize': 
        plt.figure()
        m = Basemap(llcrnrlon=-130,llcrnrlat=10,urcrnrlon=-64,urcrnrlat=50,
                projection='merc',lon_0=0)
        USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if shape_dict['HASC_1'] in corSTs.keys():
                colors[shape_dict['HASC_1']]=corSTs[shape_dict['HASC_1']]/2+0.5
            else:colors[shape_dict['HASC_1']]=0.5
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            poly = Polygon(seg,facecolor=cm.BrBG(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly)     
            
        MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                        name='states', drawbounds=True)
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if shape_dict['HASC_1'] in corSTs.keys():
                colors[shape_dict['HASC_1']]=corSTs[shape_dict['HASC_1']]/2+0.5
            else:colors[shape_dict['HASC_1']]=0.5
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            poly = Polygon(seg,facecolor=cm.BrBG(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly) 
        mapper.set_array(clrBr);plt.colorbar(mapper)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/AmENSO_yields/Yields/yieldCorr/corrMaps/'+crop+'/NorthAmerica_'+str(pvalThresh)+'.png')
        plt.close()
    elif crop is 'soy':
        plt.figure()
        m = Basemap(llcrnrlon=-140,llcrnrlat=20,urcrnrlon=-64,urcrnrlat=55,
                projection='merc',lon_0=0)
        USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if shape_dict['HASC_1'] in corSTs.keys():
                colors[shape_dict['HASC_1']]=corSTs[shape_dict['HASC_1']]/2+0.5
            else:colors[shape_dict['HASC_1']]=0.5
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            poly = Polygon(seg,facecolor=cm.BrBG(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly)
        mapper.set_array(clrBr);plt.colorbar(mapper)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/AmENSO_yields/Yields/yieldCorr/corrMaps/'+crop+'/NorthAmerica_'+str(pvalThresh)+'.png')
        plt.close() 
    elif crop is 'wheat': 
        plt.figure()
        m = Basemap(llcrnrlon=-140,llcrnrlat=10,urcrnrlon=-64,urcrnrlat=70,
                projection='merc',lon_0=0)
        USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if shape_dict['HASC_1'] in corSTs.keys():
                colors[shape_dict['HASC_1']]=corSTs[shape_dict['HASC_1']]/2+0.5
            else:colors[shape_dict['HASC_1']]=0.5
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            poly = Polygon(seg,facecolor=cm.BrBG(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly)    
     
        CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if shape_dict['HASC_1'] in corSTs.keys():
                colors[shape_dict['HASC_1']]=corSTs[shape_dict['HASC_1']]/2+0.5
            else:colors[shape_dict['HASC_1']]=0.5
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            poly = Polygon(seg,facecolor=cm.BrBG(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly)
        mapper.set_array(clrBr);plt.colorbar(mapper)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/AmENSO_yields/Yields/yieldCorr/corrMaps/'+crop+'/NorthAmerica_'+str(pvalThresh)+'.png')
        plt.close()    
       
