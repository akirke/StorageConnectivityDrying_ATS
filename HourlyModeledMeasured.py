import sys,os
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
from matplotlib import gridspec
from scipy.optimize import curve_fit
import pytz
import hydroeval as he


matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['xtick.major.pad']='12'
matplotlib.rcParams['ytick.major.pad']='8'



#Check that meteorology file looks good
import h5py
with h5py.File('shambley_met_hourly.h5', 'r') as f:
    for key in f.keys():
        print(f"{key}: {f[key].shape}")



##### Import Data from the model output file
Output = pandas.read_csv('water_balance_computational_domain.csv', skiprows=408)
Output = Output.iloc[:-1]

#Grab the times, because we don't want them in seconds
times = Output["time [s]"][:]

#And turn them into nice datetimes
start = datetime.strptime('2016-01-01 00:00:00','%Y-%m-%d %H:%M:%S')
ts = datetime.timestamp(start)
reslist = []
for sec in times:
    res_date = ts + sec
    res = datetime.fromtimestamp(res_date)
    reslist.append(res)

#Put the datetimes back in the Output Dataframe  & trim to study period
Output["DateTime"] = reslist
Output = Output.set_index('DateTime')
Output=Output.loc['2017-01-01':'2018-01-01']

#Grab the flume discharge values, convert from mol to m3, and fix the sign if needed
southFlume_out = Output["south flume discharge [mol h^-1]"][:]/(-55600.)
largeFlume_out = Output["central flume discharge [mol h^-1]"][:]/(55600.)


#Grab other water balance components
rain = Output["rain precipitation [m h^-1]"][:]
et = Output["total evapotranspiration [m h^-1]"][:]

#get to m3/ and then divide by domain size
watershed_sw_out = Output["river discharge [mol h^-1]"][:]/(55600.*1132905)
watershed_gw_out = Output["net groundwater flux [mol h^-1]"][:]/(55600.*1132905)


#Sum to Get the Monthly Water Balance
watershed_gw_month = watershed_gw_out.resample('MS').sum()
watershed_sw_month = watershed_sw_out.resample('MS').sum()
et_month = et.resample('MS').sum()
rain_month = rain.resample('MS').sum()
wb = rain_month - watershed_gw_month - watershed_sw_month - et_month

##### Plot the Monthly Water Balance
fig_wb, ax = plt.subplots(1,1, figsize=(15,10))
ax.axhline(y=0, color='k', linestyle='-')
ax.plot(et_month, 'g', label="Transpiration [Meters]")
ax.plot(watershed_sw_month, 'b', label="Surface Water Outflow [Meters]")
ax.plot(watershed_gw_month, 'brown', label="Groundwater Outflow [Meters]")
ax.plot(rain_month, 'cyan', label="Precipitation [Meters]")
ax.plot(wb, 'k', ls="--", label = "Monthly Net Water Flux")
ax.legend()
fig_wb.savefig("waterBalance.tiff", dpi=300) 


##### Get the actual flume data to compare to

#Actual 2017 Data - Central Watershed Flume
large_real_H_all = pandas.read_csv('data/LargeDischargeFlume_real.csv')
large_real_H_all=large_real_H_all.set_index('Time')
large_real_H_all.index = pandas.to_datetime(large_real_H_all.index)
large_real_H_all.index = large_real_H_all.index.tz_localize(None)


#Actual 2017 Data - South Watershed Flume
south_real_H_all = pandas.read_csv('data/SouthDischargeFlume_real.csv')
south_real_H_all=south_real_H_all.set_index('Time')
south_real_H_all.index = pandas.to_datetime(south_real_H_all.index)
south_real_H_all.index = south_real_H_all.index.tz_localize(None)


#Line up real and modeled data by their indexes
south_model_H_all = pandas.DataFrame(southFlume_out)
south_H_compare_all= pandas.merge(south_real_H_all, south_model_H_all, left_index=True, right_index=True)
south_H_compare_all = south_H_compare_all.dropna()

large_model_H_all = pandas.DataFrame(largeFlume_out)
large_H_compare_all= pandas.merge(large_real_H_all, large_model_H_all, left_index=True, right_index=True)
large_H_compare_all = large_H_compare_all.dropna()



# One to One plot of H Water Level over Whole period
fig2a, ax = plt.subplots(1,2, figsize=(30,12))

ax[0].scatter(south_H_compare_all['volume_m3'], south_H_compare_all['south flume discharge [mol h^-1]'])
ax[0].set_ylabel("Modeled")
ax[0].set_xlabel("Measured")
ax[0].set_title("Hourly Flow at South Watershed Flume [m3/h]")
ax[0].set(xlim = (0,2000))
ax[0].set(ylim = (0,2000))

ax[1].scatter(large_H_compare_all['volume_m3'], large_H_compare_all['central flume discharge [mol h^-1]'])
ax[1].set_title("Hourly Flow at Central Watershed Flume [m3/h]")
ax[1].set_ylabel("Modeled")
ax[1].set_xlabel("Measured")
ax[1].set(xlim = (0,3000))
ax[1].set(ylim = (0,3000))

fig2a.savefig("one_to_one.tiff", dpi=300) 
fig2a.show()



#####Statistics to Validate Model

## Confirmation that mean benchmark is -0.41 for modified KGE
simmean = [np.mean(south_H_compare_all["volume_m3"]) for i in range(len(south_H_compare_all["volume_m3"]))]
evaluations = south_H_compare_all["volume_m3"]

mkge, r, y, beta = he.kgeprime(simmean, evaluations.values)

print('mean benchmark mKGE= ' + str(mkge))


# mKGE for Central Watershed
evaluations = large_H_compare_all["volume_m3"]
simulations = large_H_compare_all["central flume discharge [mol h^-1]"]

mkge, r, y, beta = he.kgeprime(simulations.values, evaluations.values)

print('Central mKGE= ' + str(mkge))
print('Central mKGE r =' + str(r))
print('Central mKGE y =' + str(y))
print('Central mKGE beta =' + str(beta))


# mKGE for South Watershed
evaluations = south_H_compare_all["volume_m3"]
simulations = south_H_compare_all["south flume discharge [mol h^-1]"]

mkge, r, y, beta = he.kgeprime(simulations.values, evaluations.values)

print('South mKGE= ' + str(mkge))
print('South mKGE r =' + str(r))
print('South mKGE y =' + str(y))
print('South mKGE beta =' + str(beta))


##### Time series plot of hourly measured versus modeled discharge at weirs
#I think the super high flow day in the measured time series is a fluke
fig3a, ax = plt.subplots(2,1, figsize=(25,25))
ax[0].plot(large_model_H_all, '#d95f02', label="Modeled Flow")
ax[0].plot(large_real_H_all, 'k', label="Measured Flow")
ax[0].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
ax[0].set_ylabel('Discharge [m3/h]') 
ax[0].set(ylim=(0,3000))
ax[0].set_title("Central Watershed")
ax[0].legend()

ax[1].plot(south_model_H_all, '#d95f02', label="Modeled Flow")
ax[1].plot(south_real_H_all, 'k', label="Measured Flow")
ax[1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
ax[1].set_ylabel('Discharge [m3/h]') 
ax[1].set(ylim=(0,2000))
ax[1].set_title("South Watershed")
ax[1].legend()

fig3a.savefig("measuredModeled_ts.tiff", dpi=300) 
fig3a.show()


##### Save CSVs
south_H_compare_all.to_csv('CentralWatershed_ModeledMeasured.csv', sep=',')
large_H_compare_all.to_csv('SouthWatershed_ModeledMeasured.csv', sep=',')

