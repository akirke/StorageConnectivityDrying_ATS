#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys,os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'],'tools', 'utils'))
import ats_xdmf
import colors
import numpy as np
import matplotlib
import scipy
import scipy.stats
from matplotlib import pyplot as plt
from matplotlib import colorbar
import pandas
from datetime import datetime, date, timedelta
import matplotlib.dates as mdates
import h5py
import shapely
import shapely.geometry
import geopandas as gpd
from shapely import Polygon, MultiPolygon
import pyvista as pv
import h5py
from matplotlib.gridspec import SubplotSpec

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['xtick.major.pad']='12'
matplotlib.rcParams['ytick.major.pad']='8'




# In[2]:


#set directory as the folder where the ATS .xmf files live
directory = "VisFiles"


# In[3]:


#Load Model Visualization Files & pull out ponded depth
vis = ats_xdmf.VisFile(directory, domain='surface', output_time_unit='yr')
vis.loadMesh(order=['x','y'])
pd = vis.getArray('surface-ponded_depth')


# In[4]:


#Get Times from the days in the vis file
start = datetime.strptime('2016-01-01 00:00:00','%Y-%m-%d %H:%M:%S')
ts = datetime.timestamp(start)

times = vis.times
reslist = []

for sec in times:
    res_date = ts + (sec*24*60*60)
    res = datetime.fromtimestamp(res_date)
    reslist.append(res)


# In[5]:


#Load the center coordinates for each surface grid cell
x_centroids = vis.centroids[:,0]
y_centroids = vis.centroids[:,1]

cell_coords = list(zip(x_centroids, y_centroids))


# In[6]:


# Convert centroids to Shapely points
centroid_points = [shapely.geometry.Point(i) for i in cell_coords]


# In[7]:


# Load River Corridor Polygons North Watershed
north_Wet_gdf = gpd.read_file('StreamPoly/northWet_stream.shp')
north_Wet_geom = north_Wet_gdf.geometry.iloc[0]
if isinstance(north_Wet_geom, MultiPolygon):
    print("MultiPolygon")
    northWet = north_Wet_geom.geoms[0]
else: northWet = north_Wet_geom

north_Typ2_gdf = gpd.read_file('StreamPoly/northTyp2_stream.shp')
north_Typ2_geom = north_Typ2_gdf.geometry.iloc[0]
if isinstance(north_Typ2_geom, MultiPolygon):
    print("MultiPolygon")
    northTyp2 = north_Typ2_geom.geoms[0]
else: northTyp2 = north_Typ2_geom

north_Typ1_gdf = gpd.read_file('StreamPoly/northTyp1_stream.shp')
north_Typ1_geom = north_Typ1_gdf.geometry.iloc[0]
if isinstance(north_Typ1_geom, MultiPolygon):
    print("MultiPolygon")
    northTyp1 = north_Typ1_geom.geoms[0]
else: northTyp1 = north_Typ1_geom

#Make a mask where cells in the HGF are 'True' and the others are 'False'  
northTyp2_cells = [northTyp2.contains(centroid_point) for centroid_point in centroid_points]
northTyp1_cells = [northTyp1.contains(centroid_point) for centroid_point in centroid_points]
northWet_cells = [northWet.contains(centroid_point) for centroid_point in centroid_points]


# In[8]:


# Load River Corridor Polygons South Watershed
south_Inc_gdf = gpd.read_file('StreamPoly/southInc_stream.shp')
southInc_geom = south_Inc_gdf.geometry.iloc[0]
if isinstance(southInc_geom, MultiPolygon):
    print("MultiPolygon")
    southInc = southInc_geom.geoms[0]
else: southInc = southInc_geom

south_Typ_gdf = gpd.read_file('StreamPoly/southTyp_stream.shp')
southTyp_geom = south_Typ_gdf.geometry.iloc[0]
if isinstance(southTyp_geom, MultiPolygon):
    print("MultiPolygon")
    southTyp = southTyp_geom.geoms[0]
else: southTyp = southTyp_geom

south_Wet_gdf = gpd.read_file('StreamPoly/southWet_stream.shp')
southWet_geom = south_Wet_gdf.geometry.iloc[0]
if isinstance(southWet_geom, MultiPolygon):
    print("MultiPolygon")
    southTyp = southWet_geom.geoms[0]
else: southWet = southWet_geom

#Make a mask where cells in the HGF are 'True' and the others are 'False'  
southTyp_cells = [southTyp.contains(centroid_point) for centroid_point in centroid_points]
southWet_cells = [southWet.contains(centroid_point) for centroid_point in centroid_points]
southIncised_cells = [southInc.contains(centroid_point) for centroid_point in centroid_points]


# In[9]:


central_Inc_gdf = gpd.read_file('StreamPoly/centralInc_stream.shp')
centralInc_geom = central_Inc_gdf.geometry.iloc[0]
if isinstance(centralInc_geom, MultiPolygon):
    print("MultiPolygon")
    centralInc = centralInc_geom.geoms[0]
else: centralInc = centralInc_geom

central_Typ_gdf = gpd.read_file('StreamPoly/centralTyp_stream.shp')
centralTyp_geom = central_Typ_gdf.geometry.iloc[0]
if isinstance(centralTyp_geom, MultiPolygon):
    print("MultiPolygon")
    centralTyp = centralTyp_geom.geoms[0]
else: centralTyp = centralTyp_geom

central_Wet_gdf = gpd.read_file('StreamPoly/centralWet_stream.shp')
centralWet_geom = central_Wet_gdf.geometry.iloc[0]
if isinstance(centralWet_geom, MultiPolygon):
    print("MultiPolygon")
    centralWet = centralWet_geom.geoms[0]
else: centralWet = centralWet_geom

#Make a mask where cells in the HGF are 'True' and the others are 'False'  
centralTyp_cells = [centralTyp.contains(centroid_point) for centroid_point in centroid_points]
centralWet_cells = [centralWet.contains(centroid_point) for centroid_point in centroid_points]
centralIncised_cells = [centralInc.contains(centroid_point) for centroid_point in centroid_points]


# In[10]:


#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
northTyp1_flow = []

for i,time in enumerate(vis.times):
    water_present = []
    masked_water_northTyp1 = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(northTyp1_cells)):
        if northTyp1_cells[n] == False: 
            masked_water_northTyp1.append('NA')
        else:
            if water_present[n] == True:
                masked_water_northTyp1.append(1)
            else:
                masked_water_northTyp1.append(0)

    flow_cells = masked_water_northTyp1.count(1)
    northTyp1_flow.append(flow_cells)


#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
northTyp2_flow = []

for i,time in enumerate(vis.times):
    water_present = []
    masked_water_northTyp2 = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(northTyp2_cells)):
        if northTyp2_cells[n] == False: 
            masked_water_northTyp2.append('NA')
        else:
            if water_present[n] == True:
                masked_water_northTyp2.append(1)
            else:
                masked_water_northTyp2.append(0)

    flow_cells = masked_water_northTyp2.count(1)
    northTyp2_flow.append(flow_cells)

#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
northWet_flow = []

for i,time in enumerate(vis.times):
    water_present = []
    masked_water_northWet = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(northWet_cells)):
        if northWet_cells[n] == False: 
            masked_water_northWet.append('NA')
        else:
            if water_present[n] == True:
                masked_water_northWet.append(1)
            else:
                masked_water_northWet.append(0)

    flow_cells = masked_water_northWet.count(1)        
    northWet_flow.append(flow_cells)


# In[11]:


#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
centralTyp_flow = []
for i,time in enumerate(vis.times):
    water_present = []
    masked_water_centralTyp = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(centralTyp_cells)):
        if centralTyp_cells[n] == False: 
            masked_water_centralTyp.append('NA')
        else:
            if water_present[n] == True:
                masked_water_centralTyp.append(1)
            else:
                masked_water_centralTyp.append(0)

    flow_cells = masked_water_centralTyp.count(1)
    centralTyp_flow.append(flow_cells)

#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
centralWet_flow = []
for i,time in enumerate(vis.times):
    water_present = []
    masked_water_centralWet = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(centralWet_cells)):
        if centralWet_cells[n] == False: 
            masked_water_centralWet.append('NA')
        else:
            if water_present[n] == True:
                masked_water_centralWet.append(1)
            else:
                masked_water_centralWet.append(0)

    flow_cells = masked_water_centralWet.count(1)
    centralWet_flow.append(flow_cells)

#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water

centralIncised_flow = []

for i,time in enumerate(vis.times):
    water_present = []
    masked_water_centralIncised = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(centralIncised_cells)):
        if centralIncised_cells[n] == False: 
            masked_water_centralIncised.append('NA')
        else:
            if water_present[n] == True:
                masked_water_centralIncised.append(1)
            else:
                masked_water_centralIncised.append(0)

    flow_cells = masked_water_centralIncised.count(1)
    centralIncised_flow.append(flow_cells)


# In[12]:


#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
southTyp_flow = []

for i,time in enumerate(vis.times):
    water_present = []
    masked_water_southTyp = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(southTyp_cells)):
        if southTyp_cells[n] == False: 
            masked_water_southTyp.append('NA')
        else:
            if water_present[n] == True:
                masked_water_southTyp.append(1)
            else:
                masked_water_southTyp.append(0)

    flow_cells = masked_water_southTyp.count(1)
    southTyp_flow.append(flow_cells)

#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
southWet_flow = []

for i,time in enumerate(vis.times):
    water_present = []
    masked_water_southWet = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(southWet_cells)):
        if southWet_cells[n] == False: 
            masked_water_southWet.append('NA')
        else:
            if water_present[n] == True:
                masked_water_southWet.append(1)
            else:
                masked_water_southWet.append(0)

    flow_cells = masked_water_southWet.count(1)
    southWet_flow.append(flow_cells)

#Select Cells which are in the stream corridor for this HGF and have at least 2.5 mm of water
southIncised_flow = []

for i,time in enumerate(vis.times):
    water_present = []
    masked_water_southIncised = []
    ponded_state = pd[i,:]

    for val in ponded_state:
        if val >= 0.0025:
            water_present.append(True)
        else:
            water_present.append(False)

    for n in range(len(southIncised_cells)):
        if southIncised_cells[n] == False: 
            masked_water_southIncised.append('NA')
        else:
            if water_present[n] == True:
                masked_water_southIncised.append(1)
            else:
                masked_water_southIncised.append(0)

    flow_cells = masked_water_southIncised.count(1)
    southIncised_flow.append(flow_cells)


# In[13]:


#Calculate Percent flowing north watershed

northTyp1_extent = max(northTyp1_flow)
northTyp1_percent = [x/northTyp1_extent for x in northTyp1_flow]

northTyp2_extent = max(northTyp2_flow)
northTyp2_percent = [x/northTyp2_extent for x in northTyp2_flow]

northWet_extent = max(northWet_flow)
northWet_percent = [x/northWet_extent for x in northWet_flow]


# In[14]:


#Calculate Percent flowing south watershed
southTyp_extent = max(southTyp_flow)
southTyp_percent = [x/southTyp_extent for x in southTyp_flow]

southWet_extent = max(southWet_flow)
southWet_percent = [x/southWet_extent for x in southWet_flow]

southIncised_extent = max(southIncised_flow)
southIncised_percent = [x/southIncised_extent for x in southIncised_flow]


# In[15]:


#Calculate Percent flowing central watershed
centralTyp_extent = max(centralTyp_flow)
centralTyp_percent = [x/centralTyp_extent for x in centralTyp_flow]

centralIncised_extent = max(centralIncised_flow)
centralIncised_percent = [x/centralIncised_extent for x in centralIncised_flow]

centralWet_extent = max(centralWet_flow)
centralWet_percent = [x/centralWet_extent for x in centralWet_flow]


# In[16]:


def create_subtitle(fig: plt.Figure, grid: SubplotSpec, title: str):
    "Sign sets of subplots with title"
    row = fig.add_subplot(grid)
    # the '\n' is important
    row.set_title(f'{title}\n', fontweight='semibold')
    # hide subplot
    row.set_frame_on(False)
    row.axis('off')


# In[17]:


fig, ax = plt.subplots(3,3, figsize=(25, 25))
grid = plt.GridSpec(3, 3)
create_subtitle(fig, grid[0, ::], 'North Watershed')
create_subtitle(fig, grid[1, ::], 'South Watershed')
create_subtitle(fig, grid[2, ::], 'Central Watershed')
fig.tight_layout(pad=2)
fig.set_facecolor('w')

ax[0,0].plot(reslist, northTyp2_percent, color = '#7570b3')
ax[0,0].set_title('Upstream Intact Corridor')
ax[0,0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[0,0].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[0,0].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[0,0].set(ylim= (0,1))
ax[0,0].set_ylabel("Proportion of Channel Wet")

ax[0,1].plot(reslist, northWet_percent, color = '#1b9e77')
ax[0,1].set_title('Wetland Corridor')
ax[0,1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[0,1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[0,1].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[0,1].set(ylim= (0,1))

ax[0,2].plot(reslist, northTyp1_percent, color = '#7570b3')
ax[0,2].set_title('Downstream Intact Corridor')
ax[0,2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[0,2].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[0,2].set(ylim= (0,1))
ax[0,2].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))



ax[1,0].plot(reslist, southWet_percent, color = '#1b9e77')
ax[1,0].set_title('Wetland Corridor')
ax[1,0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[1,0].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[1,0].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[1,0].set_ylabel("Proportion of Channel Wet")
ax[1,0].set(ylim= (0,1))

ax[1,1].plot(reslist, southTyp_percent, color = '#7570b3')
ax[1,1].set_title('Intact Corridor')
ax[1,1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[1,1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[1,1].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 0:00:00','%Y-%m-%d %H:%M:%S')))
ax[1,1].set(ylim= (0,1))

ax[1,2].plot(reslist, southIncised_percent, color = '#d95f02')
ax[1,2].set_title('Incised Corridor')
ax[1,2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[1,2].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[1,2].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[1,2].set(ylim= (0,1))

ax[2,0].plot(reslist, centralTyp_percent, color = '#d95f02')
ax[2,0].set_title('Upstream Incised Corridor')
ax[2,0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[2,0].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[2,0].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[2,0].set_ylabel("Proportion of Channel Wet")
ax[2,0].set(ylim= (0,1))

ax[2,1].plot(reslist, centralWet_percent, color = '#1b9e77')
ax[2,1].set_title('Wetland Corridor')
ax[2,1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[2,1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[2,1].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[2,1].set(ylim= (0,1))

ax[2,2].plot(reslist, centralIncised_percent, color = '#d95f02')
ax[2,2].set_title('Downstream Incised Corridor')
ax[2,2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax[2,2].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[2,2].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[2,2].set(ylim= (0,1))

plt.savefig("percentFlowing.tiff", dpi=300) 


# In[18]:


#Save the df as a csv
df = pandas.DataFrame({'Date':reslist, 'NorthTyp2':northTyp2_percent, 'NorthWet': northWet_percent, 'NorthTyp1':northTyp1_percent, 'SouthWet':southWet_percent, 'SouthTyp':southTyp_percent,
                       'SouthInc':southIncised_percent, 'CentralTyp':centralTyp_percent, 'CentralWet':centralWet_percent, 'CentralInc':centralIncised_percent})

df.to_csv('percent_flowing.csv', sep=',')

