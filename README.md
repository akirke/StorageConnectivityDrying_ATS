# StorageConnectivityDrying_ATS
Workflows, output observation and visualization files, and data analysis scripts associated with the use of the Advanced Terrestrial Simulator Model for the manuscript "Storage, Connectivity, and Stream Drying in Headwater Systems: Modeling Across Hydrogeomorphic Features".

To generate input files (.xml, .exo, and .h5) an example from Watershed Workflow (Coon & Shuai, 2022) was modified (ShambleyMesh_2025.ipynb) and run in Jupyter Notebook. Input data for this workflow include: MeteorologicalData.csv and writeWeatherH5.ipynb. 

Then, ATS was run 3 times, first with constant precipitation inputs to initialize the water table (Shambley-steadystate.xml), and then a 4 year 'warm spinup' (Shambley-cyclic_steadystate.xml) was run, before the simulation used in the manuscript was initiated (Shambley-transient.xml).

Observed discharge data used to validate the model are in two csv files (ObservedDischarge_CentralWatershed and ObservedDischarge_SouthWatershed)

Data visualization and analysis based on the observation files (water_balance_computational_domain.csv) and visualization files (.xmf files) were peformed in Python, with extractPercentWet.py, HourlyModeledMeasured.py, Fig3Fig4Fig5.py, and exfil_dep_rec.py being used for figures reproduced in the main text, and DepthToWater_snapshot.py used for the creation of a supplemental figure.

Geomorphic analysis of HGFs referenced in the manuscript are based on cross sections extracted from a DEM using an R script, (shambleyCrossSections.R)
