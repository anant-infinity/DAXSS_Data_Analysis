import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter
import netCDF4 as nc
from dateutil import parser
from datetime import datetime, timedelta


plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True


# Path to Level-1 File
daxss_level1_file_path = 'C:/Users/278an/Desktop/science_analysis/301_Data/DAXSS/level1/daxss_solarSXR_level1_2022-02-14-mission_v1.0.0_test.ncdf'
# Import File as a netCDF Dataset
daxsslevel1 = nc.Dataset(daxss_level1_file_path)

time_ISO = daxsslevel1['TIME_ISO'][:]
daxss_datetime_obj_array = []
for time_index in range(len(time_ISO)):
    time_ISO_String = []
    for var_index in range(0,20):
        time_str = time_ISO[time_index][var_index].decode("utf-8")
        time_ISO_String.append(time_str)
    daxss_datetime_obj_array.append(parser.parse(''.join(time_ISO_String)))

daxss_time = date2num(daxss_datetime_obj_array)
daxss_x123_slow_count = daxsslevel1['X123_SLOW_COUNT'][:]

# Path to Level-1 File
goes_file_path = 'C:/Users/278an/Desktop/science_analysis/301_Data/GOES/sci_xrsf-l2-flx1s_g16_d20220329_v2-1-0.nc'
# Import File as a netCDF Dataset
goes_data = nc.Dataset(goes_file_path)

goes_time_array = goes_data['time'][:]
base_date = datetime(2000, 1, 1, 12, 0, 0)
goes_datetime_obj_array = []
for i in range(len(goes_time_array)):
    curr_time = base_date + timedelta(seconds=goes_time_array[i])
    goes_datetime_obj_array.append(curr_time)


goes_time = date2num(goes_datetime_obj_array)
goes_flux_data = goes_data['xrsb_flux'][:]


fig, ax1 = plt.subplots()

ax1.plot_date(daxss_time,daxss_x123_slow_count, 'o--', color='red',label="DAXSS Slow Counts")
ax1.set_ylabel("DAXSS Slow Counts (Counts/sec)",color="red",fontsize=12)
ax2 = ax1.twinx()
ax2.plot_date(goes_time,goes_flux_data, 'o:',color='blue',label="GOES XRS-B Flux")
ax2.set_ylabel("GOES XRS-B Flux (W/m^2)",color="blue",fontsize=12)

ax1.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d-%H-%M'))
ax1.tick_params(rotation=45)

plt.show()