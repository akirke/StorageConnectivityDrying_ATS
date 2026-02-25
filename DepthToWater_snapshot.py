import sys,os
sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'],'tools', 'utils'))
import ats_xdmf
import colors
import numpy as np
import matplotlib
import scipy
from scipy import stats
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
import seaborn as sns
import scikit_posthocs as sp
from matplotlib.gridspec import SubplotSpec
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['xtick.major.pad']='12'
matplotlib.rcParams['ytick.major.pad']='8'

#set directory as the folder where the ATS .xmf files live
directory = "VisFiles"


#Load Model Visualization Files & pull out depth to water
vis = ats_xdmf.VisFile(directory, domain='surface', output_time_unit='yr')
vis.loadMesh(order=['x','y'])
pd = vis.getArray('surface-water_table_depth')

#Load the center coordinates for each surface grid cell
x_centroids = vis.centroids[:,0]
y_centroids = vis.centroids[:,1]
cell_coords = list(zip(x_centroids, y_centroids))

# Convert centroids to Shapely points
centroid_points = [shapely.geometry.Point(i) for i in cell_coords]

# Load HGF Polygons from North Watershed
north_Wet_gdf = gpd.read_file('StreamPoly/northwetland.shp')
north_Wet_geom = north_Wet_gdf.geometry.iloc[0]
if isinstance(north_Wet_geom, MultiPolygon):
    print("MultiPolygon")
    northWet = north_Wet_geom.geoms[0]
else: northWet = north_Wet_geom

north_Typ2_gdf = gpd.read_file('StreamPoly/northtypical2.shp')
north_Typ2_geom = north_Typ2_gdf.geometry.iloc[0]
if isinstance(north_Typ2_geom, MultiPolygon):
    print("MultiPolygon")
    northTyp2 = north_Typ2_geom.geoms[0]
else: northTyp2 = north_Typ2_geom

north_Typ1_gdf = gpd.read_file('StreamPoly/northtypical1.shp')
north_Typ1_geom = north_Typ1_gdf.geometry.iloc[0]
if isinstance(north_Typ1_geom, MultiPolygon):
    print("MultiPolygon")
    northTyp1 = north_Typ1_geom.geoms[0]
else: northTyp1 = north_Typ1_geom


#Make a mask where cells in the HGF are 'True' and the others are 'False'  
northTyp2_cells = [northTyp2.contains(centroid_point) for centroid_point in centroid_points]
northTyp1_cells = [northTyp1.contains(centroid_point) for centroid_point in centroid_points]
northWet_cells = [northWet.contains(centroid_point) for centroid_point in centroid_points]


# Load HGF Polygons from  South Watershed
south_Inc_gdf = gpd.read_file('StreamPoly/southincised.shp')
southInc_geom = south_Inc_gdf.geometry.iloc[0]
if isinstance(southInc_geom, MultiPolygon):
    print("MultiPolygon")
    southInc = southInc_geom.geoms[0]
else: southInc = southInc_geom

south_Typ_gdf = gpd.read_file('StreamPoly/southtypical.shp')
southTyp_geom = south_Typ_gdf.geometry.iloc[0]
if isinstance(southTyp_geom, MultiPolygon):
    print("MultiPolygon")
    southTyp = southTyp_geom.geoms[0]
else: southTyp = southTyp_geom

south_Wet_gdf = gpd.read_file('StreamPoly/southwetland.shp')
southWet_geom = south_Wet_gdf.geometry.iloc[0]
if isinstance(southWet_geom, MultiPolygon):
    print("MultiPolygon")
    southTyp = southWet_geom.geoms[0]
else: southWet = southWet_geom

#Make a mask where cells in the HGF are 'True' and the others are 'False'  
southTyp_cells = [southTyp.contains(centroid_point) for centroid_point in centroid_points]
southWet_cells = [southWet.contains(centroid_point) for centroid_point in centroid_points]
southIncised_cells = [southInc.contains(centroid_point) for centroid_point in centroid_points]


# Load HGF Polygons from Central Watershed
central_Inc_gdf = gpd.read_file('StreamPoly/centralincised1.shp')
centralInc_geom = central_Inc_gdf.geometry.iloc[0]
if isinstance(centralInc_geom, MultiPolygon):
    print("MultiPolygon")
    centralInc = centralInc_geom.geoms[0]
else: centralInc = centralInc_geom

central_Typ_gdf = gpd.read_file('StreamPoly/centraltypical.shp')
centralTyp_geom = central_Typ_gdf.geometry.iloc[0]
if isinstance(centralTyp_geom, MultiPolygon):
    print("MultiPolygon")
    centralTyp = centralTyp_geom.geoms[0]
else: centralTyp = centralTyp_geom

central_Wet_gdf = gpd.read_file('StreamPoly/centralwetland.shp')
centralWet_geom = central_Wet_gdf.geometry.iloc[0]
if isinstance(centralWet_geom, MultiPolygon):
    print("MultiPolygon")
    centralWet = centralWet_geom.geoms[0]
else: centralWet = centralWet_geom

#Make a mask where cells in the HGF are 'True' and the others are 'False'  
centralTyp_cells = [centralTyp.contains(centroid_point) for centroid_point in centroid_points]
centralWet_cells = [centralWet.contains(centroid_point) for centroid_point in centroid_points]
centralIncised_cells = [centralInc.contains(centroid_point) for centroid_point in centroid_points]


#Choose the times we're interested in! Let's choose the midpoint of each season: May 5 (day 490) & Nov 5 (day 674)
timelist = vis.times.tolist()
dry_time_index = timelist.index(674)
wet_time_index = timelist.index(490)

#Pull the depth to water df for that day
dry_waterDepth = pd[dry_time_index,:]
wet_waterDepth = pd[wet_time_index,:]


#Make a list of polygon names so I don't have to type it every time
polygon_names = ["northTyp1", "northTyp2", "northWet", 
    "southIncised", "southWet", "southTyp", 
    "centralIncised", "centralTyp", "centralWet"]

#Stick all the wet/dry day masks in a results dictionary
results = {}

for name in polygon_names:
    cell_mask = globals()[f"{name}_cells"]

    #  filter depths by cell mask
    wet_list = [wet_waterDepth[n] if cell_mask[n] else np.nan for n in range(len(cell_mask))]
    dry_list = [dry_waterDepth[n] if cell_mask[n] else np.nan for n in range(len(cell_mask))]

    results[f"wet_{name}_depth"] = pandas.DataFrame(pandas.to_numeric(wet_list, errors='coerce'))
    results[f"dry_{name}_depth"] = pandas.DataFrame(pandas.to_numeric(dry_list, errors='coerce'))


#Color code and Label Code each DF column

box_pal = {"Wetland": '#1b9e77',
"Intact Riparian": '#7570b3',
    "Incised": '#d95f02'}

hgf_mapping = {
    "Wet": "Wetland",
    "Typ": "Intact Riparian",
    "Incised": "Incised"}

name_cleaner = {"northTyp2": "North Watershed Intact Riparian HGF Upstream",
    "northWet": "North Watershed Wetland HGF",
    "northTyp1": "North Watershed Intact Riparian HGF Downstream",
    "centralTyp": "Central Watershed Incised HGF Upstream",
    "centralWet": "Central Watershed Wetland HGF",
    "centralIncised": "Central Watershed Incised HGF Downstream",
    "southWet": "South Watershed Wetland HGF",
    "southTyp": "South Watershed Intact Riparian HGF",
    "southIncised": "South Watershed Incised HGF"}

#make sure they plot in the same order as the labeling
plot_order = ["North Watershed Intact Riparian HGF Upstream",
    "North Watershed Wetland HGF",
    "North Watershed Intact Riparian HGF Downstream",
    "Central Watershed Incised HGF Upstream",
    "Central Watershed Wetland HGF",
    "Central Watershed Incised HGF Downstream",
    "South Watershed Wetland HGF",
    "South Watershed Intact Riparian HGF",
    "South Watershed Incised HGF"]

#Convert Dataframe to long form
all_data = [] 
for key, df in results.items():
    #Split to get Season, Feature Name, and Width (stream vs Full HGF)
    parts = key.split('_')
    season = "Wet Season" if parts[0] == "wet" else "Dry Season"
    feature_id = parts[1] # e.g., "northTyp1"
    width = "Stream Only" if "stream" in key else "Full HGF"

    #Categorize by HGF -- remember that we redefined the upstream central watershed HGF from intact riparian to incised
    hgf_cat = "Unknown"
    if "Wet" in feature_id: hgf_cat = "Wetland"
    elif "centralTyp" in feature_id: hgf_cat = "Incised"
    elif "Typ" in feature_id: hgf_cat = "Intact Riparian"
    elif "Incised" in feature_id: hgf_cat = "Incised"

    temp_df = df.copy()
    temp_df.rename(columns={temp_df.columns[0]: 'depth'}, inplace=True)
    temp_df = temp_df.dropna(subset=['depth'])

    # Add metadata columns
    temp_df['Season'] = season
    temp_df['HGF Type'] = hgf_cat
    temp_df['Width'] = width
    temp_df['Location'] = name_cleaner.get(feature_id, feature_id) # Use clean name
    all_data.append(temp_df)

df_final = pandas.concat(all_data).dropna()

#flip the sign so water table plots below '0'
df_final['depth'] = df_final['depth'] * -1


#Filter by Season and by HGF width (full HGF versus stream)
df_full_HGF = df_final[df_final['Width'] == "Full HGF"]
df_full_HGF_wet = df_full_HGF[df_full_HGF['Season'] == "Wet Season"]
df_full_HGF_dry = df_full_HGF[df_full_HGF['Season'] == "Dry Season"]


#Statistics to check for differences across HGFs in Wet Season snapshot
data_groups = [group['depth'].values for name, group in df_full_HGF_wet.groupby('Location')]
stat, p_val = stats.kruskal(*data_groups)
print(f"Kruskal-Wallis p-value for Wet Season Snapshot: {p_val}")


dunn_results = sp.posthoc_dunn(
    df_full_HGF_wet, 
    val_col='depth', 
    group_col='Location', 
    p_adjust='Holm'
)


cld = sp.compact_letter_display(dunn_results, alpha=0.05)
print(" ")
print("Compact Letter Display (letters are different when p<0.05")
print(" ")
print(cld)


#Statistics to check for differences across HGFs in Dry Season snapshot
data_groups = [group['depth'].values for name, group in df_full_HGF_dry.groupby('Location')]
stat, p_val = stats.kruskal(*data_groups)
print(f"Kruskal-Wallis p-value for Dry Season Snapshot: {p_val}")


dunn_results = sp.posthoc_dunn(
    df_full_HGF_dry, 
    val_col='depth', 
    group_col='Location', 
    p_adjust='Holm'
)


cld = sp.compact_letter_display(dunn_results, alpha=0.05)
print(" ")
print("Compact Letter Display (letters are different when p<0.05")
print(" ")
print(cld)

#Make a figure
sns.set_style("whitegrid")
plt.figure(figsize=(14, 8))

#Pull the letters in from the Compact Letter Display which printed above
cld_map = {
    "Central Watershed Incised HGF Downstream": "a",
    "Central Watershed Incised HGF Upstream": "a",
    "Central Watershed Wetland HGF": "b",
    "North Watershed Intact Riparian HGF Downstream": "ab",
    "North Watershed Intact Riparian HGF Upstream": "ab",
    "North Watershed Wetland HGF": "b",
    "South Watershed Incised HGF": "a",
    "South Watershed Intact Riparian HGF": "a",
    "South Watershed Wetland HGF": "b"
}

ax = sns.boxplot(
    data=df_full_HGF_wet,
    x="Location",        
    y="depth",           
    hue="HGF Type",      
    palette=box_pal,     
    order=plot_order,
    showfliers=False     
)

plt.axvline(x=2.5, color='black', linestyle='--', alpha=0.3)
plt.axvline(x=5.5, color='black', linestyle='--', alpha=0.3)
plt.axhline(0, color='black', linewidth=1, alpha=0.5)

ax.set_xlabel("")      
ax.set_xticklabels([])
ax.set_ylabel("Depth to Water")      

watershed_labels = [
    (1, "NORTH WATERSHED"),
    (4, "CENTRAL WATERSHED"),
    (7, "SOUTH WATERSHED")
]


for i, location in enumerate(plot_order):
    letter = cld_map.get(location, "")

    ax.text(i, 0.005, letter, 
            ha='center', va='bottom', 
            fontsize=12, fontweight='bold', color='black')


plt.title("Water Level across HGF on May 5, 2017", fontsize=18, pad=40, fontweight='bold')
for x_coord, text in watershed_labels:
    plt.text(x_coord, 0.2, text, ha='center', va='top', 
             fontsize=12, fontweight='bold', color='#333333')

plt.savefig("May5_WL.tiff", dpi=300) 


#Make another Figure
sns.set_style("whitegrid")
plt.figure(figsize=(14, 8))

#Pull the letters in from the Compact Letter Display which printed above
cld_map = {
    "Central Watershed Incised HGF Downstream": "a",
    "Central Watershed Incised HGF Upstream": "a",
    "Central Watershed Wetland HGF": "bc",
    "North Watershed Intact Riparian HGF Downstream": "abc",
    "North Watershed Intact Riparian HGF Upstream": "abc",
    "North Watershed Wetland HGF": "c",
    "South Watershed Incised HGF": "ab",
    "South Watershed Intact Riparian HGF": "ab",
    "South Watershed Wetland HGF": "c"
}

ax = sns.boxplot(
    data=df_full_HGF_dry,
    x="Location",         
    y="depth",           
    hue="HGF Type",      
    palette=box_pal,     
    order=plot_order,
    showfliers=False     
)

plt.axvline(x=2.5, color='black', linestyle='--', alpha=0.3)
plt.axvline(x=5.5, color='black', linestyle='--', alpha=0.3)
plt.axhline(0, color='black', linewidth=1, alpha=0.5)

ax.set_xlabel("")      
ax.set_xticklabels([])
ax.set_ylabel("Depth to Water")      

watershed_labels = [
    (1, "NORTH WATERSHED"),
    (4, "CENTRAL WATERSHED"),
    (7, "SOUTH WATERSHED")
]

for i, location in enumerate(plot_order):
    letter = cld_map.get(location, "")

    ax.text(i, 0.005, letter, 
            ha='center', va='bottom', 
            fontsize=12, fontweight='bold', color='black')


plt.title("Water Level across HGF on November 5, 2017", fontsize=18, pad=40, fontweight='bold')
for x_coord, text in watershed_labels:
    plt.text(x_coord, 0.2, text, ha='center', va='top', 
             fontsize=12, fontweight='bold', color='#333333')

plt.savefig("Nov5_WL.tiff", dpi=300) 


