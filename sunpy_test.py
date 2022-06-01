import numpy as np
import sunpy.data.sample
import sunpy.timeseries as ts
import sunpy.map
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a

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
results = Fido.search(a.Time("2022/02/1", "2022/06/1"),a.Instrument.goes )
print(results.show("Start Time", 'End Time', 'Instrument', 'Physobs'))

#Downloading the raw files
#downloaded_files = Fido.fetch(results, path='C:/Users/278an/Desktop/science_analysis/301_Data/level1/{file}')