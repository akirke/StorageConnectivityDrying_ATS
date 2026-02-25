# imports data from ATS model output and plots flow from individual Hydrogeomorphic features across a watershed
# Ashleigh Kirker modified 2/24/2026

#Import packages for plotting and data analysis
import sys,os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import scipy
import scipy.stats
import pandas
from datetime import datetime, date, timedelta
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
from matplotlib.gridspec import SubplotSpec
import matplotlib.gridspec as gridspec
from sklearn.linear_model import LinearRegression
import seaborn as sns

#We are using these Color Blind Safe Category Colors
#1b9e77 - Wetland
#d95f02 - Incised
#7570b3 - Typical

directory = ""

#This is the length of streams in the watershed in meters
aims_streams = 2130.17

#Given meters of stream in each HGF, what percent of total streams does each HGF contain?
northWet_stream = 199.03
northWet_percentLength = northWet_stream/aims_streams
northTyp1_stream = 181.56
northTyp1_percentLength = northTyp1_stream/aims_streams
northTyp2_stream = 80.34
northTyp2_percentLength = northTyp2_stream/aims_streams

southTyp_stream = 157.77
southTyp_percentLength = southTyp_stream/aims_streams

southInc_stream = 352.80
southInc_percentLength = southInc_stream/aims_streams
southWet_stream = 199.53
southWet_percentLength = southWet_stream/aims_streams

centralInc_stream = 479.20
centralInc_percentLength = centralInc_stream/aims_streams
print(centralInc_percentLength)
centralIncised2_stream = 137.38
centralIncised2_percentLength = centralIncised2_stream/aims_streams
centralWet_stream = 340.71
centralWet_percentLength = centralWet_stream/aims_streams


#Load  Watershed Data from the .csv file output by ATS
Output = pandas.read_csv('water_balance_computational_domain.csv', skiprows=408).iloc[:-1] 

#Convert time from seconds to datetime
times = Output["time [s]"][:]
start = datetime.strptime('2016-01-01 00:00:00','%Y-%m-%d %H:%M:%S')
ts = datetime.timestamp(start)

reslist = []
for sec in times:
    res_date = ts + sec
    res = datetime.fromtimestamp(res_date)
    reslist.append(res)

#Set datetime as index and filter to study period
Output["DateTime"] = reslist
Output = Output.set_index('DateTime')
Output = Output.loc['2017-01-01':'2018-01-01']

#Save the filtered output to a nicer csv
Output.to_csv('2017_ModelOutput.csv', sep=',')

#Get Water Table Depths for each HGF
southIncised_depth = -Output["South Incised Reach Depth to Water [m]"][:]
southTyp_depth = -Output["South Typical Reach Depth to Water [m]"][:]
southWet_depth = -Output["South Wet Reach Depth to Water [m] "][:]

#Get Outflow from each HGF and convert from mol/h to m3/hr. Fix the Sign if needed
southIncised_out = Output["South Incised Reach Outflow [mol/h]"][:]/(55600.)
southTyp_out = Output["South Typical Reach Outflow [mol/h]"][:]/(55600.)
southWet_out = Output["South Wet Reach Outflow [mol/h]"][:]/(55600.)

#Get Water Table Depths for each HGF
northTyp1_depth = -Output["North Typical1 Reach Depth to Water [m]"][:]
northTyp2_depth = -Output["North Typical2 Reach Depth to Water [m]"][:]
northWet_depth = -Output["North Wet Reach Depth to Water [m] "][:]

#Get Outflow from each HGF and convert from mol/h to m3/hr. Fix the Sign if needed
northTyp1_out = Output["North Typical1 Reach Outflow [mol/h]"][:]/(55600.)
northTyp2_out = Output["North Typical2 Reach Outflow [mol/h]"][:]/(55600.)
northWet_out = Output["North Wet Reach Outflow [mol/h]"][:]/(55600.)


#Get Water Table Depths for each HGF
centralIncised1_depth = -Output["Central Incised Reach Depth to Water [m]"][:]
centralIncised2_depth = -Output["Central Typical Reach Depth to Water [m]"][:]
centralWet_depth = -Output["Central Wet Reach Depth to Water [m] "][:]

#Get Outflow from each HGF and convert from mol/h to m3/hr. Fix the Sign if needed
centralInc_out = Output["Central Incised Reach Outflow [mol/h]"][:]/(55600.)
centralIncised2_out = Output["Central Typical Reach Outflow [mol/h]"][:]/(55600.)
centralWet_out = Output["Central Wet Reach Outflow [mol/h]"][:]/(55600.)


###### Convert Outflows to Daily Sums. This simplifies the plots and also makes it easier to compare outflows
#from upstream and downstream HGFs without dealing with time lag

#Turn the HGF's outflow from a series into a pandas dataframe so we can manipulate it
southIncised_out_daily = pandas.DataFrame(southIncised_out)
#Sum HGF outflow across an entire day, so that we are reporting outflow in cubic meters per day
southIncised_out_daily = southIncised_out_daily.resample('D').sum()
southIncised_out_daily = southIncised_out_daily.rename(columns={0:"Runoff [m3]"})

southTyp_out_daily = pandas.DataFrame(southTyp_out)
southTyp_out_daily = southTyp_out_daily.resample('D').sum()
southTyp_out_daily = southTyp_out_daily.rename(columns={0:"Runoff [m3]"})

southWet_out_daily = pandas.DataFrame(southWet_out)
southWet_out_daily = southWet_out_daily.resample('D').sum()
southWet_out_daily = southWet_out_daily.rename(columns={"South Wet Reach Outflow [mol/h]":"Runoff [m3]"})

centralIncised1_out_daily = pandas.DataFrame(centralInc_out)
centralIncised1_out_daily = centralIncised1_out_daily.resample('D').sum()
centralIncised1_out_daily = centralIncised1_out_daily.rename(columns={0:"Runoff [m3]"})

centralIncised2_out_daily = pandas.DataFrame(centralIncised2_out)
centralIncised2_out_daily = centralIncised2_out_daily.resample('D').sum()
centralIncised2_out_daily = centralIncised2_out_daily.rename(columns={"Central Typical Reach Outflow [mol/h]":"Runoff [m3]"})

centralWet_out_daily = pandas.DataFrame(centralWet_out)
centralWet_out_daily = centralWet_out_daily.resample('D').sum()
centralWet_out_daily = centralWet_out_daily.rename(columns={0:"Runoff [m3]"})

northTyp1_out_daily = pandas.DataFrame(northTyp1_out)
northTyp1_out_daily = northTyp1_out_daily.resample('D').sum()
northTyp1_out_daily = northTyp1_out_daily.rename(columns={0:"Runoff [m3]"})

northTyp2_out_daily = pandas.DataFrame(northTyp2_out)
northTyp2_out_daily = northTyp2_out_daily.resample('D').sum()
northTyp2_out_daily = northTyp2_out_daily.rename(columns={"North Typical2 Reach Outflow [mol/h]":"Runoff [m3]"})

northWet_out_daily = pandas.DataFrame(northWet_out)
northWet_out_daily = northWet_out_daily.resample('D').sum()
northWet_out_daily = northWet_out_daily.rename(columns={0:"Runoff [m3]"})


#Instead of reporting total outflow at each HGF, it is more useful to report the flow *generated* in each HGF
#which we can calculate by subtracting the daily flow into HGF from the daily flow out of it. 
#For headwater HGFs where the stream originates, inflow is zero

northTyp2_flowgen_daily = northTyp2_out_daily.iloc[:,0]
northWet_flowgen_daily = northWet_out_daily.iloc[:,0] - northTyp2_out_daily.iloc[:,0]
northTyp1_flowgen_daily = northTyp1_out_daily.iloc[:,0] - northWet_out_daily.iloc[:,0]

centralIncised2_flowgen_daily = centralIncised2_out_daily.iloc[:,0]
centralWet_flowgen_daily = centralWet_out_daily.iloc[:,0] - centralIncised2_out_daily.iloc[:,0]
centralIncised1_flowgen_daily = centralIncised1_out_daily.iloc[:,0] - centralWet_out_daily.iloc[:,0]

southWet_flowgen_daily = southWet_out_daily.iloc[:,0]
southIncised_flowgen_daily = southIncised_out_daily.iloc[:,0] - southTyp_out_daily.iloc[:,0]
southTyp_flowgen_daily = southTyp_out_daily.iloc[:,0] - southWet_out_daily.iloc[:,0]

#The total flow generated in all of our subwatersheds is the sum of the flow leaving each subwatershed
tot_flow_daily =  centralIncised1_out_daily.iloc[:,0] + southIncised_out_daily.iloc[:,0] +  northTyp1_out_daily.iloc[:,0]


#Save all the flowgen values (negatives included) to a csv file so that we can generate donut plots in excel
flowgendf = pandas.concat([northTyp2_flowgen_daily,northWet_flowgen_daily,northTyp1_flowgen_daily, southWet_flowgen_daily, southTyp_flowgen_daily, southIncised_flowgen_daily, centralIncised2_flowgen_daily, centralWet_flowgen_daily, centralIncised1_flowgen_daily],axis=1)
flowgendf.to_csv('daily_flowgen_v2.csv', sep=',')
flowgendf.describe()

#This figure will plot flow generated in each HGF for each of three subwatersheds

#Figure 3.  Plot showing daily proportions of flow at outlet.
fig, ax = plt.subplots(3,1, figsize=(12, 20))

#stacked area plots are weird, they show negative values as positive, so we need to replace negative values with 0
northTyp2_flow_daily = northTyp2_flowgen_daily.clip(lower=0)
northWet_flow_daily = northWet_flowgen_daily.clip(lower=0)
northTyp1_flow_daily = northTyp1_flowgen_daily.clip(lower=0)
southWet_flow_daily = southWet_flowgen_daily.clip(lower=0)
southTyp_flow_daily = southTyp_flowgen_daily.clip(lower=0)
southIncised_flow_daily = southIncised_flowgen_daily.clip(lower=0)
centralIncised2_flow_daily = centralIncised2_flowgen_daily.clip(lower=0)
centralWet_flow_daily = centralWet_flowgen_daily.clip(lower=0)
centralIncised1_flow_daily = centralIncised1_flowgen_daily.clip(lower=0)

#Plot the North Watershed 
ax[0].stackplot(northTyp2_flow_daily.index,  northTyp2_flow_daily, northWet_flow_daily, northTyp1_flow_daily, labels = ['Intact Riparian HGF (upstream)', 'Wetland HGF', 'Intact Riparian HGF (downstream)'], 
                 colors=['#9770b3','#1b9e77', '#7570b3'], linewidth=0.001)
ax[0].set_title("A) Streamflow Gain from HGFs in North Watershed", fontsize=20, fontweight="bold")
ax[0].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
ax[0].legend(loc="upper right", fontsize=14)
ax[0].tick_params(axis='both', which='major', labelsize=14)



#Plot the Central Watershed 
ax[1].stackplot(northTyp2_flow_daily.index, centralIncised2_flow_daily, centralWet_flow_daily, centralIncised1_flow_daily, 
                 labels = ['Incised HGF (upstream)', 'Wetland HGF', 'Incised HGF (downstream)'], colors=['#feb37a','#1b9e77','#d95f02'], linewidth=0.001)
ax[1].set_title("C) Streamflow Gain from HGFs in Central Watershed", fontsize = 20, fontweight="bold")
ax[1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
ax[1].set_ylabel('Streamflow Gain ($m^3$ per day)', fontsize=25, labelpad=30) 
ax[1].legend(loc="upper right", fontsize=14)
ax[1].tick_params(axis='both', which='major', labelsize=14)


#Plot the South Watershed 

ax[2].stackplot(southWet_flow_daily.index, southWet_flow_daily, southTyp_flow_daily,southIncised_flow_daily, labels = ['Wetland HGF','Intact Riparain HGF', 'Incised HGF'],
              colors=['#1b9e77','#7570b3', '#d95f02'], linewidth=0.001)
ax[2].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[2].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
ax[2].set_title("B) Streamflow Gain from HGFs in South Watershed", fontsize=20, fontweight="bold")
ax[2].legend(loc="upper right", fontsize=14)
ax[2].tick_params(axis='both', which='major', labelsize=14)

#Add some space in between subplots
plt.subplots_adjust(wspace=2)
plt.savefig("Figure3.tiff", dpi=300) 


#The next portion of the code will generate a Cumulative Distribution Function Plot which aggregates data from each category of HGF
#across all subwatersheds and compares it to total outflow

#sum up generated flow from HGFs of the same category
wetland_outflow = northWet_flowgen_daily + southWet_flowgen_daily + centralWet_flowgen_daily
incised_outflow = southIncised_flowgen_daily + centralIncised1_flowgen_daily + centralIncised2_flowgen_daily
intact_outflow = northTyp2_flowgen_daily + northTyp1_flowgen_daily + southTyp_flowgen_daily 


#We are also interested in how the Stream Length in each HGF type correlates to the amount of Flow generated in each HGF
#So we will import the stream length in each HGF and calculate that HGF's proportion of total streams in the watershed

#Aggregate by HGF type
allInc_streamLength_per = centralInc_percentLength + southInc_percentLength + centralIncised2_percentLength
allWet_streamLength_per = centralWet_percentLength + southWet_percentLength + northWet_percentLength
allTyp_streamLength_per =  southTyp_percentLength + northTyp1_percentLength + northTyp2_percentLength

print(allInc_streamLength_per)
print(allWet_streamLength_per)
print(allTyp_streamLength_per)


# Consolidate summed timeseries into one long-format DataFrame

data = pandas.DataFrame({
    'Discharge': pandas.concat([incised_outflow, intact_outflow, wetland_outflow, tot_flow_daily]),
    'Category': (['Incised'] * len(incised_outflow) + 
                 ['Intact'] * len(intact_outflow) + 
                 ['Wetland'] * len(wetland_outflow) + 
                 ['Total'] * len(tot_flow_daily))
})

data_pct = pandas.DataFrame({
    'Incised Streamflow Gain_Pct': (incised_outflow / tot_flow_daily) * 100,
    'Wetland Streamflow Gain_Pct': (wetland_outflow / tot_flow_daily) * 100,
    'Intact Riparian Streamflow Gain_Pct': (intact_outflow / tot_flow_daily) * 100
})


box_pal = {
        "Wetland Streamflow Gain_Pct": '#1b9e77', "Intact Riparian Streamflow Gain_Pct" : '#7570b3',
    "Incised Streamflow Gain_Pct": '#d95f02',
}


### Plot the proportion of streamflow from each HGF at a given streamflow exceedence###
plt.figure(figsize=(12,8))

for col in ['Incised Streamflow Gain_Pct', 'Wetland Streamflow Gain_Pct', 'Intact Riparian Streamflow Gain_Pct']:
    # sort by catchment-wide outflow
    sorted_indices = tot_flow_daily.argsort()
    y_values = data_pct[col].iloc[sorted_indices]
    x_values = np.linspace(0, 1, len(y_values))
    plt.plot(x_values, y_values, label=col.replace('_Pct', ''), color=box_pal[col], lw=2)

plt.axhline(y=45, color='#d95f02', linestyle='dotted')
plt.axhline(y=35, color='#1b9e77', linestyle='dotted')
plt.axhline(y=20, color='#7570b3', linestyle='dotted')

plt.annotate(
    '% of Stream Length in Incised HGFs', #annotate to mark 45% of stream length in incised HGFs
    xy=(0, 45), 
    xytext=(10, 5), 
    textcoords='offset points',
    va='center',
    color='#642c01'
)

plt.annotate(
    '% of Stream Length in Wetland HGFs', #annotate to mark 35% of stream length in wetland HGFs
    xy=(0, 35), 
    xytext=(-30, 5), #
    textcoords='offset points',
    va='center', 
    color='#0d4a38'
)
plt.annotate(
    '% of Stream Length in Intact Riparian HGFs', #annotate to mark 20% of stream length in intact HGFs
    xy=(0, 20),
    xytext=(90, 5), 
    textcoords='offset points',
    va='center', 
    color='#2a284b'
)


plt.gca().xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter())

plt.ylabel('% Contribution to Total Discharge', fontsize = 15, fontweight="bold")
plt.xlabel('% Flow Exceedence (Dry $\\rightarrow$ Wet)', fontsize = 15, fontweight="bold")
plt.legend()
plt.show()

plt.savefig("Figure4.tiff", dpi=300) 




#Load data from the percent_flowing csv. This was generated by examining model grid cells within each HGF's stream and classifying them as 'wet' or 'dry'
#We use this %wet  metric to give an idea of the extent of the stream network at 2-day timesteps

percentFlow = pandas.read_csv('percent_flowing.csv')
percentFlow['Date'] = pandas.to_datetime(percentFlow['Date'])
percentFlow.index = percentFlow['Date']

#Aggregate data by HGF type
percentTyp = percentFlow[['NorthTyp2', 'NorthTyp1','SouthTyp']]
percentWet = percentFlow[['NorthWet', 'SouthWet', 'CentralWet']]
percentInc = percentFlow[['SouthInc', 'CentralInc','CentralTyp']]

#Calculate the average percent wet for each time step for each HGF type
avgperTyp = percentTyp.mean(axis=1)
avgperWet = percentWet.mean(axis=1)
avgperInc = percentInc.mean(axis=1)


#Plot percent flow and catchment discharge for Figure 5
fig, ax = plt.subplots(2,1, figsize=(30, 15))


#In the top panel, plot total outflow from all three subwatersheds and delineate the dry season and the wet season in beige and blue
ax[0].plot(tot_flow_daily, 'k', linewidth=4, label="Total Catchment Outfow")
ax[0].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[0].set_ylabel("Daily Discharge [$m^3$]", fontsize=30, fontweight='bold')
ax[0].axvspan(datetime.strptime('2017-03-20 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2017-06-21 00:00:00','%Y-%m-%d %H:%M:%S'), alpha=0.25, color='#01665e', label = "Wet Season")
ax[0].axvspan(datetime.strptime('2017-09-22 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2017-12-21 00:00:00','%Y-%m-%d %H:%M:%S'), alpha=0.25, color='#d8b365', label = "Dry Season")
ax[0].get_xaxis().set_visible(False)
ax[0].set(ylim=(30,30000))
ax[0].tick_params(axis='y', which='major', labelsize=30)
ax[0].legend(loc='upper right', fontsize=30)
ax[0].set_yscale('log')


#In the bottom panel, plot percent wet at each time step across watersheds, by HGF
ax[1].plot(avgperTyp,color='#7570b3', label="Intact Riparian HGFs", linewidth = 5, ls='--' )
ax[1].plot(avgperWet,color='#1b9e77', label ="Wetland HGFs", linewidth = 5, ls='--' )
ax[1].plot(avgperInc, color='#d95f02', label = "Incised HGFs",linewidth = 5, ls='--' )
ax[1].tick_params(axis='both', which='major', labelsize=30)
ax[1].set(xlim= (datetime.strptime('2017-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), datetime.strptime('2018-01-01 00:00:00','%Y-%m-%d %H:%M:%S')))
ax[1].set_ylabel("Percent of Stream Inundated", fontsize=30, fontweight='bold')
ax[1].yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
ax[1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
ax[1].legend(loc='lower right', fontsize=30)
ax[1].set(ylim=(0,1))

fig.tight_layout(h_pad=2)
plt.savefig("Figure5.tiff", dpi=300) 