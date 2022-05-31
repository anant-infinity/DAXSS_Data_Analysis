# Python Program to view Level-1 DAXSS netCDF files
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
# Path to Level-1 File
daxss_level1_file_path = 'C:/Users/278an/Desktop/science_analysis/301_Data/level1/daxss_solarSXR_level1_2022-02-14-mission_v1.0.0.ncdf'
#Import File as a netCDF Dataset
daxsslevel1 = nc.Dataset(daxss_level1_file_path)

#Printing the dataset to get information of metadata, dimensions and variables
#print(daxsslevel1)

#Creating a dictionary containing the information about the dataset
metadata_info_dict = daxsslevel1.__dict__
#Access particular element in metadata of dataset
#print (metadata_info_dict['TITLE'])

#Accessing metadata of the dimensions
for dim in daxsslevel1.dimensions.values():
    print(dim)

#Accessing metadata of the variables
for var in daxsslevel1.variables.values():
    print(var)

#Slow count
x123_slow_count = daxsslevel1['DATA.X123_SLOW_COUNT'][:]
print(x123_slow_count)

spectrum_index = 2913
#fig, ax = plt.subplots(1)
# plt.plot(daxsslevel1[spectrum_index]["energy"], daxsslevel1[spectrum_index]["irradiance"], drawstyle="steps-mid")
plt.scatter(daxsslevel1['DATA.TIME_JD'][0,:], daxsslevel1['DATA.X123_SLOW_COUNT'][0, :], c ="blue",
            marker ="s",
            edgecolor ="blue",
            s = 5)  # If loaded from netCDF
#fig.autofmt_xdate()
#ax.xaxis.set_major_formatter(mdates.DateFormatter("%b - %Y"))
plt.xlabel("Time_JD")
plt.ylim([1e3, 1e6])
plt.yscale("log")
plt.ylabel("DAXSS X123_Slow_Count (cps)")
plt.suptitle("DAXSS Solar SXR Spectrum ")
plt.show()
