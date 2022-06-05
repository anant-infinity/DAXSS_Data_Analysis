# Python Program to view Level-1 DAXSS netCDF files
import netCDF4 as nc
import matplotlib.pyplot as plt
from jdcal import jd2gcal

# Path to Level-1 File
daxss_level1_file_path = 'C:/Users/278an/Desktop/science_analysis/301_Data/DAXSS/level1/daxss_solarSXR_level1_2022-02-14-mission_v1.0.0.ncdf'
# Import File as a netCDF Dataset
daxsslevel1 = nc.Dataset(daxss_level1_file_path)

# Printing the dataset to get information of metadata, dimensions and variables
#print(daxsslevel1)

# Creating a dictionary containing the information about the dataset
metadata_info_dict = daxsslevel1.__dict__
# Access particular element in metadata of dataset
#print (metadata_info_dict['TITLE'])

# Accessing metadata of the dimensions
#for dim in daxsslevel1.dimensions.values():
#    print(dim)

# Accessing metadata of the variables
var_array = []
for var in daxsslevel1.variables.values():
    var_array.append(var)

#print(var_array[0].name)

plt.figure()
# Plotting the X123 Slow Count
plt.subplot(121)
x123_slow_count = daxsslevel1['DATA.X123_SLOW_COUNT'][:]
plt.scatter(daxsslevel1['DATA.TIME_JD'][0,:], daxsslevel1['DATA.X123_SLOW_COUNT'][0, :], c ="blue",
            marker ="s",
            edgecolor ="blue",
            s = 5)

plt.xlabel("Time_JD")
plt.ylim([1e3, 1e6])
plt.yscale("log")
plt.ylabel("DAXSS X123_Slow_Count (cps)")
plt.suptitle("DAXSS X123 Slow Count ")

# Plotting Spectrum
plt.subplot(122)
index_quiet_sun = 0
irradiance_quite_sun = daxsslevel1['DATA.IRRADIANCE'][0,0,:]
index_flare = 1381
irradiance_flare = daxsslevel1['DATA.IRRADIANCE'][0,1381,:]
index_eclipse = 3245
irradiance_eclipse = daxsslevel1['DATA.IRRADIANCE'][0,3245,:]
energy = daxsslevel1['DATA.ENERGY'][0,0,:]

plt.plot(energy, irradiance_quite_sun, color='k', label='quite sun')
plt.plot(energy, irradiance_flare, color='r', label='flare')
plt.plot(energy, irradiance_eclipse, color='b', label='eclipse')

plt.xlabel("Energy (keV)")
plt.ylim([1e2, 1e9])
plt.yscale("log")
plt.ylabel("Irradiance (ph/s/cm^2/keV)")
plt.suptitle("DAXSS Example Spectra ")
plt.legend()
plt.show()

jd2gcal(2000,1,1)