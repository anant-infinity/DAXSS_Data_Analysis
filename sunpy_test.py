import numpy as np
import sunpy.data.sample
import sunpy.timeseries as ts
import sunpy.map
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
import netCDF4 as nc
# Plotting an example time series
my_timeseries = ts.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
#my_timeseries.peek()

# Plotting an example MAP
aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
fig = plt.figure()
ax = plt.subplot(111, projection=aia)
aia.plot()
aia.draw_limb()
aia.draw_grid()
plt.colorbar()
#plt.show()

# Searching for Datasets
#print(a.Instrument)

# Searching for GOES XRS Data
#results = Fido.search(a.Time("2022/05/27", "2022/06/1"),a.Instrument.goes )
#print(results.show("Start Time", 'End Time', 'Instrument', 'Physobs'))

#Downloading the raw files
#downloaded_files = Fido.fetch(results, path='C:/Users/278an/Desktop/science_analysis/301_Data/GOES/{file}')




# Path to Level-1 File
goes_file_path = 'C:/Users/278an/Desktop/science_analysis/301_Data/GOES/sci_xrsf-l2-flx1s_g16_d20220527_v2-1-0.nc'
# Import File as a netCDF Dataset
goes_data = nc.Dataset(goes_file_path)

# Printing the dataset to get information of metadata, dimensions and variables
# print(goes_data)

# Accessing metadata of the dimensions
for dim in goes_data.dimensions.values():
    print(dim)

# Accessing metadata of the variables
for var in goes_data.variables.values():
    print(var)

plt.figure()
# Plotting the X123 Slow Count

plt.plot(goes_data['time'][:], goes_data['xrsa_flux'][:])

plt.xlabel("Time_LS")
plt.ylim([1e-7, 1e-5])
plt.yscale("log")
plt.ylabel("GOES Primary XRS-A channel flux")
plt.suptitle("Primary XRS-A channel flux")
plt.show()