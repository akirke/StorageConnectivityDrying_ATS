# StorageConnectivityDrying_ATS
Workflows, output observation and visualization files, and data analysis scripts associated with the use of the Advanced Terrestrial Simulator Model for the manuscript "Storage, Connectivity, and Stream Drying in Headwater Systems: Modeling Across Hydrogeomorphic Features".

To generate input files (xml and h5) an example from Watershed Workflow (Coon & Shuai, 2022) was modified ().
Then, ATS was run 3 times, first with constant precipitation inputs to initialize the water table, and then a 4 year 'warm spinup' (cyclic_steadystate) was run, before the simulation used in the manuscript was initiated (transient).

Data visualization and analysis based on the observation files (water-balance-computational-domain.csv) and visualization files (.xmf files in VisFiles.zip) were peformed in Python, with extractPercentWet.py, HourlyModeledMeasured.py, Fig3Fig4Fig5.py, and exfil_dep_rec.py being used for figures reproduced in the main text, and DepthToWater_snapshot.py used for the creation of a supplemental figure.

Geomorphic analysis of HGFs referenced in the manuscript are based on cross sections extracted from a DEM using an R script, (shambleyCrossSections.R)
