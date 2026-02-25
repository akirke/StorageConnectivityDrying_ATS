import sys,os
import numpy as np
import matplotlib
import scipy
from scipy import stats
from matplotlib import pyplot as plt
import pandas
from datetime import datetime, date, timedelta
import matplotlib.dates as mdates
import h5py
from matplotlib import gridspec
import shapely
import shapely.geometry
from scipy.optimize import curve_fit
import scikit_posthocs as sp
import piecewise_regression
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick

#HGU sizes in square meters
southIncised_size = 14371.71 # 45%
southTyp_size = 5674.57
southWet_size = 11700.50
all_south_HGF = 31746.78 # 

centralIncised1_size = 16044.52 # 44%
centralWet_size = 15386.22
centralIncised2_size = 5110.05 #upstream
all_central_HGF = 36540.79

northTyp1_size = 8188.58
northTyp2_size = 2193.07
northWet_size = 7254.86
all_north_HGF = 18149.51

southWat_size = 118456.04
centralWat_size = 264960.39
northWat_size = 110165.0


#Load output from big observations csv file 
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


##### South Watershed HGFs #####
#Convert Outflow mol/h to m3/h
southIncised_out = Output["South Incised Reach Outflow [mol/h]"][:]/(55600.)
southTyp_out = Output["South Typical Reach Outflow [mol/h]"][:]/(55600.)
southWet_out = Output["South Wet Reach Outflow [mol/h]"][:]/(55600.)

#Flip the sign so depth to water plots better
southIncised_depth = -Output["South Incised Reach Depth to Water [m]"][:]
southTyp_depth = -Output["South Typical Reach Depth to Water [m]"][:]
southWet_depth = -Output["South Wet Reach Depth to Water [m] "][:]

#Convert Exfiltration mol/h to m3/h
southIncised_exfil = Output["South Incised Reach Exfiltration [mol/h]"][:]/(55600.)
southTyp_exfil = Output["South Typical Reach Exfiltration [mol/h]"][:]/(55600.)
southWet_exfil = Output["South Wet Reach Exfiltration [mol/h]"][:]/(55600.)


##### Central Watershed HGFs #####
#Convert Outflow mol/h to m3/h
centralIncised_out = Output["Central Incised Reach Outflow [mol/h]"][:]/(55600.)
centralIncised2_out = Output["Central Typical Reach Outflow [mol/h]"][:]/(55600.) #This one has been reclassified from intact to incised
centralWet_out = Output["Central Wet Reach Outflow [mol/h]"][:]/(55600.)

#Flip the sign so depth to water plots better
centralIncised_depth = -Output["Central Incised Reach Depth to Water [m]"][:]
centralIncised2_depth = -Output["Central Typical Reach Depth to Water [m]"][:]
centralWet_depth = -Output["Central Wet Reach Depth to Water [m] "][:]

#Convert Exfiltration mol/h to m3/h
centralIncised_exfil = Output["Central Incised Reach Exfiltration [mol/h]"][:]/(55600.)
centralIncised2_exfil = Output["Central Typical Reach Exfiltration [mol/h]"][:]/(55600.)
centralWet_exfil = Output["Central Wet Reach Exfiltration [mol/h]"][:]/(55600.)


##### North Watershed HGFs #####
#Convert Outflow mol/h to m3/h
northTyp1_out = Output["North Typical1 Reach Outflow [mol/h]"][:]/(55600.)
northTyp2_out = Output["North Typical2 Reach Outflow [mol/h]"][:]/(55600.)
northWet_out = Output["North Wet Reach Outflow [mol/h]"][:]/(55600.)

#Flip the sign so depth to water plots better
northTyp1_depth = -Output["North Typical1 Reach Depth to Water [m]"][:]
northTyp2_depth = -Output["North Typical2 Reach Depth to Water [m]"][:]
northWet_depth = -Output["North Wet Reach Depth to Water [m] "][:]

#Convert Exfiltration mol/h to m3/h
northTyp1_exfil = Output["North Typical1 Reach Exfiltration [mol/h]"][:]/(55600.)
northTyp2_exfil = Output["North Typical2 Reach Exfiltration [mol/h]"][:]/(55600.)
northWet_exfil = Output["North Wet Reach Exfiltration [mol/h]"][:]/(55600.)



#Get change in storage in millimeters instead of meters
southIncised_exfil_mm = (southIncised_exfil/southIncised_size)*1000
southTyp_exfil_mm = (southTyp_exfil/southTyp_size)*1000
southWet_exfil_mm = (southWet_exfil/southWet_size)*1000

centralIncised_exfil_mm = (centralIncised_exfil/centralIncised1_size)*1000
centralIncised2_exfil_mm = (centralIncised2_exfil/centralIncised2_size)*1000
centralWet_exfil_mm = (centralWet_exfil/centralWet_size)*1000

northTyp1_exfil_mm = (northTyp1_exfil/northTyp1_size)*1000
northTyp2_exfil_mm = (northTyp2_exfil/northTyp2_size)*1000
northWet_exfil_mm = (northWet_exfil/northWet_size)*1000

#Calculate Daily Change in Storage (gain/loss) for each reach, then get go from Daily to a median weekly value to reduce autocorrelation
southIncised_exfil_daily = pandas.DataFrame(southIncised_exfil_mm)
southIncised_exfil_daily = southIncised_exfil_daily.resample('D').sum()
southIncised_exfil_daily = southIncised_exfil_daily.resample('W').median()
southIncised_exfil_daily = southIncised_exfil_daily.rename(columns={'South Incised Reach Exfiltration [mol/h]':'a'})

southTyp_exfil_daily = pandas.DataFrame(southTyp_exfil_mm)
southTyp_exfil_daily = southTyp_exfil_daily.resample('D').sum()
southTyp_exfil_daily = southTyp_exfil_daily.resample('W').median()
southTyp_exfil_daily = southTyp_exfil_daily.rename(columns={'South Typical Reach Exfiltration [mol/h]':'a'})

southWet_exfil_daily = pandas.DataFrame(southWet_exfil_mm)
southWet_exfil_daily = southWet_exfil_daily.resample('D').sum()
southWet_exfil_daily = southWet_exfil_daily.resample('W').median()
southWet_exfil_daily = southWet_exfil_daily.rename(columns={'South Wet Reach Exfiltration [mol/h]':'a'})

centralIncised_exfil_daily = pandas.DataFrame(centralIncised_exfil_mm)
centralIncised_exfil_daily = centralIncised_exfil_daily.resample('D').sum()
centralIncised_exfil_daily = centralIncised_exfil_daily.resample('W').median()
centralIncised_exfil_daily = centralIncised_exfil_daily.rename(columns={'Central Incised Reach Exfiltration [mol/h]':'a'})

centralIncised2_exfil_daily = pandas.DataFrame(centralIncised2_exfil_mm)
centralIncised2_exfil_daily = centralIncised2_exfil_daily.resample('D').sum()
centralIncised2_exfil_daily = centralIncised2_exfil_daily.resample('W').median()
centralIncised2_exfil_daily = centralIncised2_exfil_daily.rename(columns={'Central Typical Reach Exfiltration [mol/h]':'a'})

centralWet_exfil_daily = pandas.DataFrame(centralWet_exfil_mm)
centralWet_exfil_daily = centralWet_exfil_daily.resample('D').sum()
centralWet_exfil_daily = centralWet_exfil_daily.resample('W').median()
centralWet_exfil_daily = centralWet_exfil_daily.rename(columns={'Central Wet Reach Exfiltration [mol/h]':'a'})

northTyp1_exfil_daily = pandas.DataFrame(northTyp1_exfil_mm)
northTyp1_exfil_daily = northTyp1_exfil_daily.resample('D').sum()
northTyp1_exfil_daily = northTyp1_exfil_daily.resample('W').median()
northTyp1_exfil_daily = northTyp1_exfil_daily.rename(columns={'North Typical1 Reach Exfiltration [mol/h]':'a'})

northTyp2_exfil_daily = pandas.DataFrame(northTyp2_exfil_mm)
northTyp2_exfil_daily = northTyp2_exfil_daily.resample('D').sum()
northTyp2_exfil_daily = northTyp2_exfil_daily.resample('W').median()
northTyp2_exfil_daily = northTyp2_exfil_daily.rename(columns={'North Typical2 Reach Exfiltration [mol/h]':'a'})

northWet_exfil_daily = pandas.DataFrame(northWet_exfil_mm)
northWet_exfil_daily = northWet_exfil_daily.resample('D').sum()
northWet_exfil_daily = northWet_exfil_daily.resample('W').median()
northWet_exfil_daily = northWet_exfil_daily.rename(columns={'North Wet Reach Exfiltration [mol/h]':'a'})


#Change from Outflow in m3/hr to average daily m3/hr
southIncised_out_daily = pandas.DataFrame(southIncised_out)
southIncised_out_daily = southIncised_out_daily.resample('D').mean()
southTyp_out_daily = pandas.DataFrame(southTyp_out)
southTyp_out_daily = southTyp_out_daily.resample('D').mean()
southWet_out_daily = pandas.DataFrame(southWet_out)
southWet_out_daily = southWet_out_daily.resample('D').mean()

centralIncised_out_daily = pandas.DataFrame(centralIncised_out)
centralIncised_out_daily = centralIncised_out_daily.resample('D').mean()
centralIncised2_out_daily = pandas.DataFrame(centralIncised2_out)
centralIncised2_out_daily = centralIncised2_out_daily.resample('D').mean()
centralWet_out_daily = pandas.DataFrame(centralWet_out)
centralWet_out_daily = centralWet_out_daily.resample('D').mean()

northTyp1_out_daily = pandas.DataFrame(northTyp1_out)
northTyp1_out_daily = northTyp1_out_daily.resample('D').mean()
northTyp2_out_daily = pandas.DataFrame(northTyp2_out)
northTyp2_out_daily = northTyp2_out_daily.resample('D').mean()
northWet_out_daily = pandas.DataFrame(northWet_out)
northWet_out_daily = northWet_out_daily.resample('D').mean()

#Calculate streamflow gain per HGF by subtracting outflow minus inflow
north_Q_wet_daily = northWet_out_daily.iloc[:,0] - northTyp2_out_daily.iloc[:,0]
north_Q_wet_daily = pandas.DataFrame(north_Q_wet_daily)
north_Q_wet_daily = north_Q_wet_daily.rename(columns={0:"Q"})
north_Q_typ2_daily = northTyp2_out_daily.iloc[:,0]
north_Q_typ2_daily = pandas.DataFrame(north_Q_typ2_daily)
north_Q_typ2_daily = north_Q_typ2_daily.rename(columns={"North Typical2 Reach Outflow [mol/h]":"Q"})
north_Q_typ1_daily = northTyp1_out_daily.iloc[:,0] - northWet_out_daily.iloc[:,0]
north_Q_typ1_daily = pandas.DataFrame(north_Q_typ1_daily)
north_Q_typ1_daily = north_Q_typ1_daily.rename(columns={0:"Q"})

central_Q_wet_daily = centralWet_out_daily.iloc[:,0] - centralIncised2_out_daily.iloc[:,0]
central_Q_wet_daily = pandas.DataFrame(central_Q_wet_daily)
central_Q_wet_daily = central_Q_wet_daily.rename(columns={0:"Q"})
central_Q_inc2_daily = centralIncised2_out_daily.iloc[:,0]
central_Q_inc2_daily = pandas.DataFrame(central_Q_inc2_daily)
central_Q_inc2_daily = central_Q_inc2_daily.rename(columns={"Central Typical Reach Outflow [mol/h]":"Q"})
central_Q_inc_daily = centralIncised_out_daily.iloc[:,0] - centralWet_out_daily.iloc[:,0]
central_Q_inc_daily = pandas.DataFrame(central_Q_inc_daily)
central_Q_inc_daily = central_Q_inc_daily.rename(columns={0:"Q"})


south_Q_inc_daily = southIncised_out_daily.iloc[:,0] - southTyp_out_daily.iloc[:,0]
south_Q_inc_daily = pandas.DataFrame(south_Q_inc_daily)
south_Q_inc_daily = south_Q_inc_daily.rename(columns={0:"Q"})
south_Q_wet_daily = southWet_out_daily.iloc[:,0]
south_Q_wet_daily = pandas.DataFrame(south_Q_wet_daily)
south_Q_wet_daily = south_Q_wet_daily.rename(columns={"South Wet Reach Outflow [mol/h]":"Q"})
south_Q_typ_daily = southTyp_out_daily.iloc[:,0] - southWet_out_daily.iloc[:,0]
south_Q_typ_daily = pandas.DataFrame(south_Q_typ_daily)
south_Q_typ_daily = south_Q_typ_daily.rename(columns={0:"Q"})

########### START OF RECESSION MATH ############
#Calculate daily dQ/dT
#dT is always 1 day

south_Q_wet_daily["dQ"] = south_Q_wet_daily.diff()
south_Q_typ_daily["dQ"] = south_Q_typ_daily.diff()
south_Q_inc_daily["dQ"] = south_Q_inc_daily.diff()

central_Q_wet_daily["dQ"] = central_Q_wet_daily.diff()
central_Q_inc2_daily["dQ"] = central_Q_inc2_daily.diff()
central_Q_inc_daily["dQ"] = central_Q_inc_daily.diff()

north_Q_typ1_daily["dQ"] = north_Q_typ1_daily.diff()
north_Q_typ2_daily["dQ"] = north_Q_typ2_daily.diff()
north_Q_wet_daily["dQ"] = north_Q_wet_daily.diff()

#Calculate d(dQ) --> where this value is negative, dQ is decreasing, and where dQ is negative, Q is decreasing
south_Q_wet_daily["ddQ"]= south_Q_wet_daily["dQ"].diff()
south_Q_typ_daily["ddQ"]= south_Q_typ_daily["dQ"].diff()
south_Q_inc_daily["ddQ"]= south_Q_inc_daily["dQ"].diff()

central_Q_wet_daily["ddQ"] = central_Q_wet_daily["dQ"].diff()
central_Q_inc2_daily["ddQ"] = central_Q_inc2_daily["dQ"].diff()
central_Q_inc_daily["ddQ"] = central_Q_inc_daily["dQ"].diff()

north_Q_typ1_daily["ddQ"] = north_Q_typ1_daily["dQ"].diff()
north_Q_typ2_daily["ddQ"] = north_Q_typ2_daily["dQ"].diff()
north_Q_wet_daily["ddQ"] = north_Q_wet_daily["dQ"].diff()


#Select days where both dQ and Q are decreasing 
south_wet_decreasing = south_Q_wet_daily[south_Q_wet_daily['ddQ'] > 0]
south_wet_decreasing = south_wet_decreasing[south_wet_decreasing['dQ'] < 0]
south_wet_decreasing = south_wet_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
south_wet_decreasing = south_wet_decreasing[south_wet_decreasing['Q'] > 0.0]

south_typ_decreasing = south_Q_typ_daily[south_Q_typ_daily['ddQ'] < 0]
south_typ_decreasing = south_typ_decreasing[south_typ_decreasing['dQ'] < 0]
south_typ_decreasing = south_typ_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
south_typ_decreasing = south_typ_decreasing[south_typ_decreasing['Q'] > 0.0]

south_inc_decreasing = south_Q_inc_daily[south_Q_inc_daily['ddQ'] > 0]
south_inc_decreasing = south_inc_decreasing[south_inc_decreasing['dQ'] < 0]
south_inc_decreasing = south_inc_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
south_inc_decreasing = south_inc_decreasing[south_inc_decreasing['Q'] > 0.0]

central_wet_decreasing = central_Q_wet_daily[central_Q_wet_daily['ddQ'] > 0]
central_wet_decreasing = central_wet_decreasing[central_wet_decreasing['dQ'] < 0]
central_wet_decreasing = central_wet_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
central_wet_decreasing = central_wet_decreasing[central_wet_decreasing['Q'] > 0.0]

central_inc2_decreasing = central_Q_inc2_daily[central_Q_inc2_daily['ddQ'] > 0]
central_inc2_decreasing = central_inc2_decreasing[central_inc2_decreasing['dQ'] < 0]
central_inc2_decreasing = central_inc2_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
central_inc2_decreasing = central_inc2_decreasing[central_inc2_decreasing['Q'] > 0.0]

central_inc_decreasing = central_Q_inc_daily[central_Q_inc_daily['ddQ'] > 0]
central_inc_decreasing = central_inc_decreasing[central_inc_decreasing['dQ'] < 0]
central_inc_decreasing = central_inc_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
central_inc_decreasing = central_inc_decreasing[central_inc_decreasing['Q'] > 0.0]

north_typ1_decreasing = north_Q_typ1_daily[north_Q_typ1_daily['ddQ'] > 0]
north_typ1_decreasing = north_typ1_decreasing[north_typ1_decreasing['dQ'] < 0]
north_typ1_decreasing = north_typ1_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
north_typ1_decreasing = north_typ1_decreasing[north_typ1_decreasing['Q'] > 0.0]

north_typ2_decreasing = north_Q_typ2_daily[north_Q_typ2_daily['ddQ'] > 0]
north_typ2_decreasing = north_typ2_decreasing[north_typ2_decreasing['dQ'] < 0]
north_typ2_decreasing = north_typ2_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
north_typ2_decreasing = north_typ2_decreasing[north_typ2_decreasing['Q'] > 0.0]

north_wet_decreasing = north_Q_wet_daily[north_Q_wet_daily['ddQ'] > 0]
north_wet_decreasing = north_wet_decreasing[north_wet_decreasing['dQ'] < 0]
north_wet_decreasing = north_wet_decreasing.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
north_wet_decreasing = north_wet_decreasing[north_wet_decreasing['Q'] > 0.0]



#Log Transform Q and dQ values
log_central_inc_Q = np.log(central_inc_decreasing['Q'])
log_central_inc_dQ = np.log(-1*central_inc_decreasing['dQ'])

log_central_wet_Q = np.log(central_wet_decreasing['Q'])
log_central_wet_dQ = np.log(-1*central_wet_decreasing['dQ'])

log_central_inc2_Q = np.log(central_inc2_decreasing['Q'])
log_central_inc2_dQ = np.log(-1*central_inc2_decreasing['dQ'])

log_south_inc_Q = np.log(south_inc_decreasing['Q'])
log_south_inc_dQ = np.log(-1*south_inc_decreasing['dQ'])

log_south_wet_Q = np.log(south_wet_decreasing['Q'])
log_south_wet_dQ = np.log(-1*south_wet_decreasing['dQ'])

log_south_typ_Q = np.log(south_typ_decreasing['Q'])
log_south_typ_dQ = np.log(-1*south_typ_decreasing['dQ'])

log_north_wet_Q = np.log(north_wet_decreasing['Q'])
log_north_wet_dQ = np.log(-1*north_wet_decreasing['dQ'])

log_north_typ1_Q = np.log(north_typ1_decreasing['Q'])
log_north_typ1_dQ = np.log(-1*north_typ1_decreasing['dQ'])

log_north_typ2_Q = np.log(north_typ2_decreasing['Q'])
log_north_typ2_dQ = np.log(-1*north_typ2_decreasing['dQ'])

#Group recession df by HGF typ
log_typ_Q = pandas.concat([log_north_typ1_Q, log_south_typ_Q], ignore_index=True)
log_typ_Q = pandas.concat([log_typ_Q, log_north_typ2_Q], ignore_index=True)     
log_typ_dQ = pandas.concat([log_north_typ1_dQ, log_south_typ_dQ], ignore_index=True)
log_typ_dQ = pandas.concat([log_typ_dQ, log_north_typ2_dQ], ignore_index=True)
typ_frames = [log_typ_Q, log_typ_dQ]

log_wet_Q = pandas.concat([log_north_wet_Q, log_central_wet_Q], ignore_index=True)
log_wet_Q = pandas.concat([log_wet_Q, log_south_wet_Q], ignore_index=True)
log_wet_dQ = pandas.concat([log_north_wet_dQ, log_central_wet_dQ], ignore_index=True)
log_wet_dQ = pandas.concat([log_wet_dQ, log_south_wet_dQ], ignore_index=True)
wet_frames = [log_wet_Q, log_wet_dQ]

log_inc_Q = pandas.concat([log_central_inc_Q, log_south_inc_Q], ignore_index=True)
log_inc_Q = pandas.concat([log_inc_Q, log_central_inc2_Q], ignore_index=True)
log_inc_dQ = pandas.concat([log_central_inc_dQ, log_south_inc_dQ], ignore_index=True)
log_inc_dQ = pandas.concat([log_inc_dQ, log_central_inc2_dQ], ignore_index=True)
inc_frames = [log_inc_Q, log_inc_dQ]

log_inc_Q = pandas.concat(inc_frames, axis=1)
log_wet_Q = pandas.concat(wet_frames, axis=1)
log_typ_Q = pandas.concat(typ_frames, axis=1)

log_inc_Q = log_inc_Q.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
log_wet_Q = log_wet_Q.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
log_typ_Q = log_typ_Q.replace([np.inf, -np.inf], np.nan).dropna(axis=0)

#Calcualate recession relationship slopes 
inc_b, inc_a, inc_r, inc_p, inc_err  = scipy.stats.linregress(log_inc_Q["Q"], log_inc_Q["dQ"])
wet_b, wet_a, wet_r, wet_p, wet_err  = scipy.stats.linregress(log_wet_Q["Q"], log_wet_Q["dQ"])
typ_b, typ_a, typ_r, typ_p, typ_err  = scipy.stats.linregress(log_typ_Q["Q"], log_typ_Q["dQ"])

xseq = np.linspace(-7,7, num=100)


#Check whether the slopes significantly different
def compare_slopes(b1, se1, b2, se2):
    # Calculate Z-score
    z = (b1 - b2) / np.sqrt(se1**2 + se2**2)
    # Calculate two-tailed p-value
    p_val = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p_val

z_stat, p_val = compare_slopes(wet_b, wet_err,
                               inc_b, inc_err)

print(f"Wet versus Incised Recession Slope Z-score: {z_stat:.4f}, p-value: {p_val:.5f}")

z_stat, p_val = compare_slopes(wet_b, wet_err, 
                               typ_b, typ_err)
print(f"Wet versus Intact Recession Slope Z-score: {z_stat:.4f}, p-value: {p_val:.5f}")

z_stat, p_val = compare_slopes(typ_b, typ_err, 
                               inc_b, inc_err)
print(f"Intact versus Incised Recession Slope Z-score: {z_stat:.4f}, p-value: {p_val:.10f}")

######## START OF EXFILTRATION CALCULATION ########
##### Get exfiltration values for just the wet season, spring equinox to summer solstice #####
winter_southIncised_exfil_daily = southIncised_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_southIncised_exfil_daily.reset_index(drop=True, inplace=True)
winter_southWet_exfil_daily = southWet_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_southWet_exfil_daily.reset_index(drop=True, inplace=True)
winter_southTyp_exfil_daily = southTyp_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_southTyp_exfil_daily.reset_index(drop=True, inplace=True)

winter_centralIncised_exfil_daily = centralIncised_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_centralIncised_exfil_daily.reset_index(drop=True, inplace=True)
winter_centralIncised2_exfil_daily = centralIncised2_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_centralIncised2_exfil_daily.reset_index(drop=True, inplace=True)
winter_centralWet_exfil_daily = centralWet_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_centralWet_exfil_daily.reset_index(drop=True, inplace=True)

winter_northWet_exfil_daily = northWet_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_northWet_exfil_daily.reset_index(drop=True, inplace=True)
winter_northTyp1_exfil_daily = northTyp1_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_northTyp1_exfil_daily.reset_index(drop=True, inplace=True)
winter_northTyp2_exfil_daily = northTyp2_exfil_daily.loc['2017-03-20':'2017-06-21']
winter_northTyp2_exfil_daily.reset_index(drop=True, inplace=True)

##### Get exfiltration values for just the dry season, fall equinox to winter solstice #####
summer_southIncised_exfil_daily = southIncised_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_southIncised_exfil_daily.reset_index(drop=True, inplace=True)
summer_southWet_exfil_daily = southWet_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_southWet_exfil_daily.reset_index(drop=True, inplace=True)
summer_southTyp_exfil_daily = southTyp_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_southTyp_exfil_daily.reset_index(drop=True, inplace=True)

summer_centralIncised_exfil_daily = centralIncised_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_centralIncised_exfil_daily.reset_index(drop=True, inplace=True)
summer_centralIncised2_exfil_daily = centralIncised2_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_centralIncised2_exfil_daily.reset_index(drop=True, inplace=True)
summer_centralWet_exfil_daily = centralWet_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_centralWet_exfil_daily.reset_index(drop=True, inplace=True)

summer_northWet_exfil_daily = northWet_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_northWet_exfil_daily.reset_index(drop=True, inplace=True)
summer_northTyp1_exfil_daily = northTyp1_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_northTyp1_exfil_daily.reset_index(drop=True, inplace=True)
summer_northTyp2_exfil_daily = northTyp2_exfil_daily.loc['2017-09-22':'2017-12-21']
summer_northTyp2_exfil_daily.reset_index(drop=True, inplace=True)

##### Save daily exfiltration dfs to csvs #####
southIncised_exfil_daily.to_csv('southInc_exfil.csv', sep=',')
southWet_exfil_daily.to_csv('southWet_exfil.csv', sep=',')
southTyp_exfil_daily.to_csv('southTyp_exfil.csv', sep=',')

centralIncised_exfil_daily.to_csv('centralInc1_exfil.csv', sep=',')
centralIncised2_exfil_daily.to_csv('centralInc2_exfil.csv', sep=',')
centralWet_exfil_daily.to_csv('centralWet_exfil.csv', sep=',')

northWet_exfil_daily.to_csv('northWet_exfil.csv', sep=',')
northTyp1_exfil_daily.to_csv('northTyp1_exfil.csv', sep=',')
northTyp2_exfil_daily.to_csv('northTyp2_exfil.csv', sep=',')

##### Group exfiltration dfs by season and HGF type #####
winterTyp = pandas.concat([winter_northTyp1_exfil_daily, winter_northTyp2_exfil_daily], ignore_index=True)
winterTyp = pandas.concat([winterTyp, winter_southTyp_exfil_daily], ignore_index=True)

winterWet = pandas.concat([winter_northWet_exfil_daily, winter_centralWet_exfil_daily], ignore_index=True)
winterWet = pandas.concat([winterWet, winter_southWet_exfil_daily], ignore_index=True)

winterIncised = pandas.concat([winter_centralIncised_exfil_daily, winter_southIncised_exfil_daily],ignore_index=True)
winterIncised = pandas.concat([winterIncised, winter_centralIncised2_exfil_daily],ignore_index=True)

summerTyp = pandas.concat([summer_northTyp1_exfil_daily, summer_southTyp_exfil_daily], ignore_index=True)
summerTyp = pandas.concat([summerTyp, summer_northTyp2_exfil_daily], ignore_index=True)

summerWet = pandas.concat([summer_northWet_exfil_daily, summer_centralWet_exfil_daily], ignore_index=True)
summerWet = pandas.concat([summerWet, summer_southWet_exfil_daily], ignore_index=True)

summerIncised = pandas.concat([summer_centralIncised_exfil_daily, summer_southIncised_exfil_daily] ,ignore_index=True)
summerIncised = pandas.concat([summerIncised, summer_centralIncised2_exfil_daily] ,ignore_index=True)


#Get mean and stdev exfiltration for each hgf and season
seasonframes = [ winterWet, summerWet, winterTyp, summerTyp, winterIncised, summerIncised]
seasonnames = [ "Wet Season Wetland Corridor", "Dry Season Wetland Corridor", "Wet Season Intact Riparian Corridor", "Dry Season Intact Riparian Corridor", "Wet Season Incised Corridor","Dry Season Incised Corridor"]
all_season_df = pandas.concat(seasonframes, axis=1, ignore_index=True)
all_season_df.columns = seasonnames
print(all_season_df.describe())

#Check to see whether exfiltration values are significantly different
stat, p_val = stats.kruskal(*seasonframes)
print(f"Kruskal-Wallis p-value for Exfiltration: {p_val}")

order=["Wet Season Wetland Corridor", "Dry Season Wetland Corridor", "Wet Season Intact Riparian Corridor", 
    "Dry Season Intact Riparian Corridor", "Wet Season Incised Corridor","Dry Season Incised Corridor"]

df_long = all_season_df.melt(
    value_vars=order,
    var_name='Group', 
    value_name='Exfiltration'
)

dunn_results = sp.posthoc_dunn(
    df_long, 
    val_col='Exfiltration', 
    group_col='Group', 
    p_adjust='holm'
)

dunn_sorted = dunn_results.loc[order, order]
print(" ")
print(dunn_sorted)
print(" ")
cld = sp.compact_letter_display(dunn_sorted, alpha=0.05)
cld_ordered = cld.reindex(order)

print(cld_ordered)

####### START OF WATER TABLE DEPTH CALCULATIONS #######
#Convert Water Table Depth to Daily Measurements
reslist = Output.index

#Calculate take the weekly median of water table depth so that we reduce autocorrelation
southIncised_depth_weekly = pandas.DataFrame(southIncised_depth)
southIncised_depth_weekly = southIncised_depth_weekly.resample('W').median()
southIncised_depth_weekly = southIncised_depth_weekly.rename(columns={'South Incised Reach Depth to Water [m]':'a'})

southTyp_depth_weekly = pandas.DataFrame(southTyp_depth)
southTyp_depth_weekly = southTyp_depth_weekly.resample('W').median()
southTyp_depth_weekly = southTyp_depth_weekly.rename(columns={'South Typical Reach Depth to Water [m]':'a'})

southWet_depth_weekly = pandas.DataFrame(southWet_depth)
southWet_depth_weekly = southWet_depth_weekly.resample('W').median()
southWet_depth_weekly = southWet_depth_weekly.rename(columns={'South Wet Reach Depth to Water [m] ':'a'})

centralIncised_depth_weekly = pandas.DataFrame(centralIncised_depth)
centralIncised_depth_weekly = centralIncised_depth_weekly.resample('W').median()
centralIncised_depth_weekly = centralIncised_depth_weekly.rename(columns={'Central Incised Reach Depth to Water [m]':'a'})

centralIncised2_depth_weekly = pandas.DataFrame(centralIncised2_depth)
centralIncised2_depth_weekly = centralIncised2_depth_weekly.resample('W').median()
centralIncised2_depth_weekly = centralIncised2_depth_weekly.rename(columns ={'Central Typical Reach Depth to Water [m]':'a'})

centralWet_depth_weekly = pandas.DataFrame(centralWet_depth)
centralWet_depth_weekly = centralWet_depth_weekly.resample('W').median()
centralWet_depth_weekly = centralWet_depth_weekly.rename(columns={'Central Wet Reach Depth to Water [m] ':'a'})

northTyp1_depth_weekly = pandas.DataFrame(northTyp1_depth)
northTyp1_depth_weekly = northTyp1_depth_weekly.resample('W').median()
northTyp1_depth_weekly = northTyp1_depth_weekly.rename(columns={'North Typical1 Reach Depth to Water [m]':'a'})

northTyp2_depth_weekly = pandas.DataFrame(northTyp2_depth)
northTyp2_depth_weekly = northTyp2_depth_weekly.resample('W').median()
northTyp2_depth_weekly = northTyp2_depth_weekly.rename(columns={'North Typical2 Reach Depth to Water [m]':'a'})

northWet_depth_weekly = pandas.DataFrame(northWet_depth)
northWet_depth_weekly = northWet_depth_weekly.resample('W').median()
northWet_depth_weekly = northWet_depth_weekly.rename(columns={'North Wet Reach Depth to Water [m] ':'a'})

##### Get water table depth values for just the wet season, spring equinox to summer solstice #####
winter_southIncised_depth_weekly = southIncised_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_southIncised_depth_weekly.reset_index(drop=True, inplace=True)
winter_southWet_depth_weekly = southWet_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_southWet_depth_weekly.reset_index(drop=True, inplace=True)
winter_southTyp_depth_weekly = southTyp_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_southTyp_depth_weekly.reset_index(drop=True, inplace=True)

winter_centralIncised_depth_weekly = centralIncised_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_centralIncised_depth_weekly.reset_index(drop=True, inplace=True)
winter_centralIncised2_depth_weekly = centralIncised2_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_centralIncised2_depth_weekly.reset_index(drop=True, inplace=True)
winter_centralWet_depth_weekly = centralWet_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_centralWet_depth_weekly.reset_index(drop=True, inplace=True)

winter_northWet_depth_weekly = northWet_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_northWet_depth_weekly.reset_index(drop=True, inplace=True)
winter_northTyp1_depth_weekly = northTyp1_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_northTyp1_depth_weekly.reset_index(drop=True, inplace=True)
winter_northTyp2_depth_weekly = northTyp2_depth_weekly.loc['2017-03-20':'2017-06-21']
winter_northTyp2_depth_weekly.reset_index(drop=True, inplace=True)

##### Get water table depth values for just the dry season, fall equinox to winter solstice #####
summer_southIncised_depth_weekly = southIncised_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_southIncised_depth_weekly.reset_index(drop=True, inplace=True)
summer_southWet_depth_weekly = southWet_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_southWet_depth_weekly.reset_index(drop=True, inplace=True)
summer_southTyp_depth_weekly = southTyp_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_southTyp_depth_weekly.reset_index(drop=True, inplace=True)

summer_centralIncised_depth_weekly = centralIncised_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_centralIncised_depth_weekly.reset_index(drop=True, inplace=True)
summer_centralIncised2_depth_weekly = centralIncised2_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_centralIncised2_depth_weekly.reset_index(drop=True, inplace=True)
summer_centralWet_depth_weekly = centralWet_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_centralWet_depth_weekly.reset_index(drop=True, inplace=True)

summer_northWet_depth_weekly = northWet_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_northWet_depth_weekly.reset_index(drop=True, inplace=True)
summer_northTyp1_depth_weekly = northTyp1_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_northTyp1_depth_weekly.reset_index(drop=True, inplace=True)
summer_northTyp2_depth_weekly = northTyp2_depth_weekly.loc['2017-09-22':'2017-12-21']
summer_northTyp2_depth_weekly.reset_index(drop=True, inplace=True)

#Concatenate DFs by season and HGF
winterTyp_depth = pandas.concat([winter_northTyp1_depth_weekly, winter_northTyp2_depth_weekly], ignore_index=True)
winterTyp_depth = pandas.concat([winterTyp_depth, winter_southTyp_depth_weekly], ignore_index=True)

winterWet_depth = pandas.concat([winter_northWet_depth_weekly, winter_centralWet_depth_weekly], ignore_index=True)
winterWet_depth = pandas.concat([winterWet_depth, winter_southWet_depth_weekly], ignore_index=True)

winterIncised_depth = pandas.concat([winter_centralIncised_depth_weekly, winter_southIncised_depth_weekly],ignore_index=True)
winterIncised_depth = pandas.concat([winterIncised_depth, winter_centralIncised2_depth_weekly],ignore_index=True)

summerTyp_depth = pandas.concat([summer_northTyp1_depth_weekly, summer_southTyp_depth_weekly], ignore_index=True)
summerTyp_depth = pandas.concat([summerTyp_depth, summer_northTyp2_depth_weekly], ignore_index=True)

summerWet_depth = pandas.concat([summer_northWet_depth_weekly, summer_centralWet_depth_weekly], ignore_index=True)
summerWet_depth = pandas.concat([summerWet_depth, summer_southWet_depth_weekly], ignore_index=True)

summerIncised_depth = pandas.concat([summer_centralIncised_depth_weekly, summer_southIncised_depth_weekly], ignore_index=True)
summerIncised_depth = pandas.concat([summerIncised_depth, summer_centralIncised2_depth_weekly], ignore_index=True)

#Get mean and stdev by season and HGF
depthframes = [winterWet_depth, summerWet_depth, winterTyp_depth, summerTyp_depth, winterIncised_depth, summerIncised_depth]
seasonnames = ["Wet Season Wetland Corridor", "Dry Season Wetland Corridor", "Wet Season Intact Riparian Corridor", "Dry Season Intact Riparian Corridor",  "Wet Season Incised Corridor","Dry Season Incised Corridor"]
all_depth_df = pandas.concat(depthframes, axis=1, ignore_index=True)
all_depth_df.columns = seasonnames
print(all_depth_df.describe())


#Check to see whether depth is statistically different by HGF and season
stat, p_val = stats.kruskal(*depthframes)
print(f"Kruskal-Wallis p-value for Water Table Depth: {p_val}")

order = ["Wet Season Wetland Corridor", "Dry Season Wetland Corridor", "Wet Season Intact Riparian Corridor", 
    "Dry Season Intact Riparian Corridor", "Wet Season Incised Corridor","Dry Season Incised Corridor"]

df_long = all_depth_df.melt(
    value_vars=order,
    var_name='Group', 
    value_name='Depth'
)

dunn_results = sp.posthoc_dunn(
    df_long, 
    val_col='Depth', 
    group_col='Group', 
    p_adjust='holm'
)

dunn_sorted = dunn_results.loc[order, order]
print(" ")
print(dunn_sorted)
print(" ")
cld = sp.compact_letter_display(dunn_sorted, alpha=0.05)
cld_ordered = cld.reindex(order)

print(cld)


##### CREATE FIGURE 2 #####
box_pal = {
        "Wet Season Wetland Corridor": '#1b9e77', "Dry Season Wetland Corridor": 'mintcream',
"Wet Season Intact Riparian Corridor": '#7570b3', "Dry Season Intact Riparian Corridor": 'lavender',
    "Wet Season Incised Corridor": '#d95f02', "Dry Season Incised Corridor": 'mistyrose'
}

fig = plt.figure(figsize=(30, 20))
gs = gridspec.GridSpec(2, 6, height_ratios=[3, 2])

# Plot exfiltration by season and reach
ax1 = plt.subplot(gs[0, 0:3])
sns.boxplot(data=all_season_df, palette=box_pal, showfliers=False, ax=ax1)
sns.despine(top=True)
sns.set(font_scale=1.7)
sns.set_style("white")
ax1.axhline(y=0, lw=2, color='k')
ax1.set_xticklabels([])
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.set_ylim(-2,20)
ax1.set_ylabel("Daily Exfiltration Flux [mm]", fontsize=30)
#statistical significance letters
letters = ['a', 'b', 'ac', 'bd', 'c', 'd'] 
# Position them slightly above the top whisker
for i, col in enumerate(all_season_df.columns):
    q3 = all_season_df[col].quantile(0.975)
    ax1.text(i, q3 + .8, letters[i], 
             ha='center', va='bottom', fontsize=22)

ax1.set_title("Exfiltration by Season and HGF", fontsize=30, fontweight='bold')
ax1.text(-0.17, 0.95, 'A)', transform=ax1.transAxes, fontsize=35, fontweight='bold', verticalalignment='top')


# Plot depth to water table by season and reach
ax2 = plt.subplot(gs[0, 3:6]) 
sns.boxplot(data=all_depth_df, palette=box_pal, showfliers=False, ax=ax2)
sns.despine(top=True)
ax2.set_xticklabels([])
ax2.tick_params(axis='both', which='major', labelsize=25)
ax2.set_ylim(-1.2,0)
ax2.set_ylabel("Weekly Water Table Depth [m]", fontsize=30)
letters = ['a', 'b', 'b', 'c', 'b', 'c']
for i, col in enumerate(all_depth_df.columns):
    q3 = all_depth_df[col].quantile(0.975)
    ax2.text(i, q3 + 0.08, letters[i], 
             ha='center', va='bottom', fontsize=22)
ax2.set_title("Depth to Water Table by Season and HGF", fontsize=30, fontweight='bold')
ax2.text(-0.17, 0.95, 'B)', transform=ax2.transAxes, fontsize=35, fontweight='bold', verticalalignment='top')


# Recession plot showing dry down in each reach 
ax3 =  plt.subplot(gs[1,0:2])
ax3.plot(log_wet_Q["Q"], log_wet_Q["dQ"], '#1b9e77', marker='.', ls='', markersize=15, alpha=0.4)
ax3.plot(xseq, wet_a + wet_b * xseq, color='#1b9e77', lw=4, label=f"Slope = {wet_b:.2f}")
ax3.set_xlabel("log(Q) [$\\mathregular {m^3 h^{-1}}$]", fontsize=30)
ax3.set_ylabel(r'log($\frac{-dQ}{dt})$ [$\mathregular {m^3 h^{-2}}$]', fontsize=30)      
ax3.tick_params(axis='both', which='major', labelsize=25)
ax3.set_ylim(-7.9, 7.9)
ax3.set_xlim(-8, 8) 
ax3.legend(loc="upper right", fontsize=30)
ax3.set_title("Recession: Wetland HGFs", fontweight='bold',fontsize=35)
ax3.text(-0.2, 0.95, 'C-3)', transform=ax3.transAxes, fontsize=35, fontweight='bold', verticalalignment='top')

ax4=  plt.subplot(gs[1,2:4])
ax4.plot(log_typ_Q["Q"], log_typ_Q["dQ"], '#7570b3', marker='.', ls='', markersize=15, alpha=0.4)
ax4.plot(xseq, typ_a + typ_b * xseq, color='#7570b3', lw=4, label=f"Slope = {typ_b:.2f}")
ax4.set_xlabel("log(Q) [$\\mathregular {m^3 h^{-1}}$]", fontsize=30)
ax4.set_ylabel(r'log($\frac{-dQ}{dt})$ [$\mathregular {m^3 h^{-2}}$]', fontsize=30)         
ax4.tick_params(axis='both', which='major', labelsize=25)
ax4.set_ylim(-7.9, 7.9)
ax4.set_xlim(-8, 8)
ax4.legend(loc="upper right", fontsize=30)
ax4.set_title("Recession: Intact Riparian HGFs", fontweight='bold',fontsize=35)
ax4.text(-0.2, 0.95, 'C-2)', transform=ax4.transAxes, fontsize=35, fontweight='bold', verticalalignment='top')


ax5 = plt.subplot(gs[1,4:6])
ax5.plot(log_inc_Q["Q"], log_inc_Q["dQ"], '#d95f02', marker='.', ls='', markersize=15, alpha=0.4)
ax5.plot(xseq, inc_a + inc_b * xseq, color='#d95f02', lw=4, label=f"Slope = {inc_b:.2f}")
ax5.set_xlabel("log(Q) [$\\mathregular {m^3 h^{-1}}$]", fontsize=30)
ax5.set_ylabel(r'log($\frac{-dQ}{dt})$ [$\mathregular {m^3 h^{-2}}$]', fontsize=30)        
ax5.tick_params(axis='both', which='major', labelsize=25)
ax5.set_ylim(-7.9,7.9)
ax5.set_xlim(-8, 8)
ax5.legend(loc="upper right", fontsize=30)
ax5.set_title("Recession: Incised HGFs", fontweight='bold',fontsize=35)
ax5.text(-0.2, 0.95, 'C-3)', transform=ax5.transAxes, fontsize=35, fontweight='bold', verticalalignment='top')
fig.tight_layout(pad=2)

plt.savefig("Figure2.tiff", dpi=300) 


##### BONUS CALCULATIONS ######

#Check for differences across HGFs, not accounting for seasons
TypEx = pandas.concat([northTyp1_exfil_daily, southTyp_exfil_daily], ignore_index=True)
TypEx = pandas.concat([TypEx, northTyp2_exfil_daily], ignore_index=True)

WetEx = pandas.concat([northWet_exfil_daily, centralWet_exfil_daily], ignore_index=True)
WetEx = pandas.concat([WetEx, southWet_exfil_daily], ignore_index=True)

IncisedEx = pandas.concat([centralIncised_exfil_daily, southIncised_exfil_daily] ,ignore_index=True)
IncisedEx = pandas.concat([IncisedEx, centralIncised2_exfil_daily], ignore_index=True)

allex=pandas.concat([WetEx, TypEx] ,ignore_index=True)
allex=pandas.concat([allex, IncisedEx] ,ignore_index=True)

ExFrames = [  IncisedEx, WetEx, TypEx, ]
ExNames = ["Incised Corridor", "Wetland Corridor", "Intact Riparian Corridor" ]
all_ex_df = pandas.concat(ExFrames, axis=1, ignore_index=True)
all_ex_df.columns = ExNames
print(all_ex_df.describe())

stat, p = stats.kruskal(all_ex_df['Wetland Corridor'], all_ex_df['Incised Corridor'], all_ex_df['Intact Riparian Corridor'])
print(f"Exfiltration Kruskal-Wallis p-value: {p}")

df_long = all_ex_df.melt(
    value_vars=['Wetland Corridor', 'Incised Corridor', 'Intact Riparian Corridor'], 
    var_name='Group', 
    value_name='Exfiltration'
)

dunn_results = sp.posthoc_dunn(
    df_long, 
    val_col='Exfiltration', 
    group_col='Group', 
    p_adjust='holm'
)

print(dunn_results)


#Check for differences across HGFs, not accounting for seasons
TypDepth = pandas.concat([northTyp1_depth_weekly, southTyp_depth_weekly], ignore_index=True)
TypDepth = pandas.concat([TypDepth, northTyp2_depth_weekly], ignore_index=True)

WetDepth = pandas.concat([northWet_depth_weekly, centralWet_depth_weekly], ignore_index=True)
WetDepth = pandas.concat([WetDepth, southWet_depth_weekly], ignore_index=True)

IncisedDepth = pandas.concat([centralIncised_depth_weekly, southIncised_depth_weekly] ,ignore_index=True)
IncisedDepth = pandas.concat([IncisedDepth, centralIncised2_depth_weekly], ignore_index=True)

ExFrames = [  IncisedDepth, WetDepth, TypDepth ]
ExNames = ["Incised Corridor", "Wetland Corridor", "Intact Riparian Corridor" ]

all_dep_df = pandas.concat(ExFrames, axis=1, ignore_index=True)
all_dep_df.columns = ExNames
print(all_dep_df.describe())

stat, p = stats.kruskal(all_dep_df['Wetland Corridor'], all_dep_df['Incised Corridor'], all_dep_df['Intact Riparian Corridor'])
print(f"Depth Kruskal-Wallis p-value: {p}")

df_long_dep = all_dep_df.melt(
    value_vars=['Wetland Corridor', 'Incised Corridor', 'Intact Riparian Corridor'], 
    var_name='Group', 
    value_name='Depth'
)

dunn_results = sp.posthoc_dunn(
    df_long_dep, 
    val_col='Depth', 
    group_col='Group', 
    p_adjust='holm'
)

print(dunn_results)