from xspec import *
import matplotlib.pyplot as plt
import netCDF4 as nc
from astropy.io import fits
from matplotlib.dates import date2num, DateFormatter
import numpy as np
import os
import shutil

def generateAvgSpec(daxsslevel1, meas_index):
    print('Hello')
    global mean_quiet_sun, mean_quiet_sun_error

    channel_number_array = []
    # quality_array = []
    systematic_error_array = []
    for i in range(1, 1001, 1):
        channel_number_array.append(np.int32(i))
        # quality_array.append(np.int16(1) - np.int16(daxsslevel1['VALID_FLAG'][meas_index,i+5]))
        systematic_error_array.append(np.float32(
            daxsslevel1['SPECTRUM_CPS_ACCURACY'][meas_index, i + 5] / daxsslevel1['SPECTRUM_CPS'][meas_index, i + 5]))
    c1 = channel_number_array

    mean_quiet_sun = (daxsslevel1['SPECTRUM_CPS'][meas_index, 6:1006] +
                      daxsslevel1['SPECTRUM_CPS'][meas_index + 1, 6:1006] +
                      daxsslevel1['SPECTRUM_CPS'][meas_index + 2, 6:1006] +
                      daxsslevel1['SPECTRUM_CPS'][meas_index + 3, 6:1006]) / 4

    mean_quiet_sun_error = ((daxsslevel1['SPECTRUM_CPS_PRECISION'][meas_index, 6:1006] ** 2 +
                             daxsslevel1['SPECTRUM_CPS_PRECISION'][meas_index + 1, 6:1006] ** 2 +
                             daxsslevel1['SPECTRUM_CPS_PRECISION'][meas_index + 2, 6:1006] ** 2 +
                             daxsslevel1['SPECTRUM_CPS_PRECISION'][meas_index + 3, 6:1006] ** 2) ** 0.5) / 4

    c4 = systematic_error_array  # Accuracy = Systematic Error
    # c5 = quality_array
    # Creating and Storing the FITS File
    time_ISO_String = []
    for var_index in range(0, 20):
        time_str = daxsslevel1['TIME_ISO'][int(meas_index)][var_index].decode("utf-8")
        time_ISO_String.append(time_str)

    hdr_dummy['FILENAME'] = 'minxss_fm3_avg_' + ''.join(time_ISO_String).replace(':', '-') + '.pha'
    hdr_dummy['DATE'] = ''.join(time_ISO_String).replace(':', '-')

    hdr_data['FILENAME'] = hdr_dummy['FILENAME']
    hdr_data['DATE'] = hdr_dummy['DATE']

    # Data
    hdu_data = fits.BinTableHDU.from_columns(
        [fits.Column(name='CHANNEL', format='J', array=c1),
         fits.Column(name='RATE', format='E', array=mean_quiet_sun),
         fits.Column(name='STAT_ERR', format='E', array=mean_quiet_sun_error),
         fits.Column(name='SYS_ERR', format='E', array=c4)], header=hdr_data)
    # fits.Column(name='QUALITY', format='J', array=c5)],
    dummy_primary = fits.PrimaryHDU(header=hdr_dummy)
    hdul = fits.HDUList([dummy_primary, hdu_data])
    #
    filename = pha_path + '/minxss_fm3_avg_' + ''.join(time_ISO_String).replace(':', '-') + '.pha'
    hdul.writeto(filename, overwrite=True)

    avg_fit_param_array = fitModel(filename)
    return avg_fit_param_array
    #print(avg_fit_param_array)

def generateFITSFile(daxsslevel1,meas_index,subAverage):
    #Creating the FITS file for DAXSS

    channel_number_array = []
    #quality_array = []
    systematic_error_array = []
    for i in range(1,1001,1):
        channel_number_array.append(np.int32(i))
        #quality_array.append(np.int16(1) - np.int16(daxsslevel1['VALID_FLAG'][meas_index,i+5]))
        systematic_error_array.append(np.float32(daxsslevel1['SPECTRUM_CPS_ACCURACY'][meas_index, i+5]/daxsslevel1['SPECTRUM_CPS'][meas_index,i+5]))
    c1 = channel_number_array

    # Modify the spectrum to subtract the mean of first four spectra from the actual specta
    if (subAverage == 1):

        mod_spec = daxsslevel1['SPECTRUM_CPS'][meas_index,6:1006] - mean_quiet_sun
        for i in range (0,len(mod_spec),1):
            if(mod_spec[i]<0):
                mod_spec[i]=0
        c2 = mod_spec
        c3 = (daxsslevel1['SPECTRUM_CPS_PRECISION'][meas_index,6:1006]**2 + mean_quiet_sun_error**2)**0.5 # Precision = Statitical Error
    elif(subAverage ==0):
        c2 = daxsslevel1['SPECTRUM_CPS'][meas_index, 6:1006]
        c3 = daxsslevel1['SPECTRUM_CPS_PRECISION'][meas_index, 6:1006]

    c4 = systematic_error_array  # Accuracy = Systematic Error
    #c5 = quality_array
    # Creating and Storing the FITS File
    time_ISO_String = []
    for var_index in range(0, 20):
        time_str = daxsslevel1['TIME_ISO'][int(meas_index)][var_index].decode("utf-8")
        time_ISO_String.append(time_str)

    hdr_dummy['FILENAME'] = 'minxss_fm3_PHA_'+''.join(time_ISO_String).replace(':', '-')+'.pha'
    hdr_dummy['DATE'] = ''.join(time_ISO_String).replace(':', '-')

    hdr_data['FILENAME'] = hdr_dummy['FILENAME']
    hdr_data['DATE'] =  hdr_dummy['DATE']

    # Data
    hdu_data = fits.BinTableHDU.from_columns(
            [fits.Column(name='CHANNEL', format='J', array=c1),
             fits.Column(name='RATE', format='E', array=c2),
             fits.Column(name='STAT_ERR', format='E', array=c3),
             fits.Column(name='SYS_ERR', format='E', array=c4)],header=hdr_data)
             #fits.Column(name='QUALITY', format='J', array=c5)],
    dummy_primary = fits.PrimaryHDU(header=hdr_dummy)
    hdul = fits.HDUList([dummy_primary, hdu_data])
    #
    filename = pha_path+'/minxss_fm3_PHA_'+''.join(time_ISO_String).replace(':', '-')+'.pha'
    hdul.writeto(filename, overwrite=True)


    current_param_array = []
    current_param_array.append(meas_index)
    current_param_array.append(''.join(time_ISO_String))
    fit_param_array = fitModel(filename)
    current_param_array+= fit_param_array
    main_param_array.append(current_param_array)

def fitModel(filename):
    # Clearing Old Data + Models
    AllData.clear()
    AllModels.clear()

    # Setting Feldman Abundances
    Xset.abund = 'file FITS_Files/feld_extd'
    logFile = Xset.openLog(logs_path+'/'+filename[-39:-4]+'_Log.txt')
    Xset.chatter = 8
    #Xset.show()
    # Loading the DAXSS spectrum
    spec = Spectrum(filename)
    spec.ignore('**-0.7 4.0-**')

    # define the model
    m1 = Model("vvapec")

    # Free some parameters that are frozen (Mg, Si, and S)
    #m1.vvapec.Ne.frozen = False
    m1.vvapec.Mg.frozen = False
    m1.vvapec.Si.frozen = False
    m1.vvapec.S.frozen = False
    m1.vvapec.Ca.frozen = False
    m1.vvapec.Fe.frozen = False

    # do the fit
    Fit.nIterations = 100
    Fit.query = 'yes'
    Fit.perform()
    print("Statistic is:",Fit.statistic)
    print("DOF is: ",Fit.dof)
    print("Critical Ch-sq Delta is:", Fit.criticalDelta)
    print("Reduced Chi-sq is: ",Fit.statistic/Fit.dof)

    Plot.xAxis = 'keV'
    Plot('ld', 'delc')
    ene = Plot.x(plotGroup=1, plotWindow=1)
    eneErr = Plot.xErr(plotGroup=1, plotWindow=1)
    spec = Plot.y(plotGroup=1, plotWindow=1)
    specErr = Plot.yErr(plotGroup=1, plotWindow=1)

    fitmodel = Plot.model(plotGroup=1, plotWindow=1)

    delchi = Plot.y(plotGroup=1, plotWindow=2)
    delchiErr = Plot.yErr(plotGroup=1, plotWindow=2)

    # Plotting the FIT
    fig0 = plt.figure(num=None, figsize=(12.4, 9.6), facecolor='w', edgecolor='k')
    ax0 = fig0.add_axes([0.1, 0.4, 0.8, 0.55])
    ax0.xaxis.set_visible(False)
    plt.title("DAXSS Spectum ON: " + filename[-24:-4])
    plt.errorbar(ene, spec, xerr=eneErr, yerr=specErr, fmt='.', ms=0.5, capsize=1.0, lw=0.8, label='Data')
    plt.step(ene, fitmodel, where='mid', label='Model')
    plt.yscale("log")
    plt.xlim([0.7, 4.0])
    # plt.ylim([1, 1e6])
    plt.legend()
    plt.ylabel('Rate (counts s$^{-1}$ keV$^{-1}$)')

    ax1 = fig0.add_axes([0.1, 0.1, 0.8, 0.3])
    plt.axhline(0, linestyle='dashed', color='black')
    plt.errorbar(ene, delchi, xerr=eneErr, yerr=delchiErr, fmt='.', ms=0.1, capsize=1.0, lw=0.8)
    plt.xlim([0.7, 4.0])
    plt.ylabel('$\Delta \chi$ (model-data/error)')
    plt.xlabel('Energy (keV)')
    plt.savefig(plots_path+'/'+filename[-39:-4]+"_Fit_Plot.png")
    #plt.show()
    plt.close()
    #Store FIT Parameter Results
    fit_param_results = []

    array_kT = [m1.vvapec.kT.name] + m1.vvapec.kT.values
    fit_param_results.append(array_kT)

    array_Ne = [m1.vvapec.Ne.name] + m1.vvapec.Ne.values
    fit_param_results.append(array_Ne)

    array_Mg = [m1.vvapec.Mg.name] + m1.vvapec.Mg.values
    fit_param_results.append(array_Mg)

    array_Si = [m1.vvapec.Si.name] + m1.vvapec.Si.values
    fit_param_results.append(array_Si)

    array_S = [m1.vvapec.S.name] + m1.vvapec.S.values
    fit_param_results.append(array_S)

    array_Ca = [m1.vvapec.Ca.name] + m1.vvapec.Ca.values
    fit_param_results.append(array_Ca)

    array_Fe = [m1.vvapec.Fe.name] + m1.vvapec.Fe.values
    fit_param_results.append(array_Fe)

    array_norm = [m1.vvapec.norm.name] + m1.vvapec.norm.values
    fit_param_results.append(array_norm)

    array_chi_sq = ["Red. Chi Sq."] + [Fit.statistic/Fit.dof]
    fit_param_results.append(array_chi_sq)

    return fit_param_results
    
def generateParamPlot(daxsslevel1, subAverage, avg_fit_param_array):

    time_str_array = []
    T_array = []
    Ne_array = []
    Mg_array = []
    Si_array = []
    S_array = []
    Ca_array = []
    Fe_array = []
    norm_array = []
    red_chi_sq_array = []

    Error_T_array = []
    Error_Ne_array = []
    Error_Mg_array = []
    Error_Si_array = []
    Error_S_array = []
    Error_Ca_array = []
    Error_Fe_array = []
    Error_norm_array = []


    daxss_x123_slow_count = []
    for i in range (len(main_param_array)):
        #Ignoring data points with absurdly large abundances or reduced chi-squared values
        if(main_param_array[i][3][1] > 3 or main_param_array[i][4][1] > 3 or main_param_array[i][5][1] > 3 or main_param_array[i][6][1] > 3
            or main_param_array[i][7][1] > 3 or main_param_array[i][8][1] > 3 or main_param_array[i][10][1] > 5 or main_param_array[i][10][1] < 0.2 ):
            continue
        time_str_array.append(main_param_array[i][1])
        T_array.append(main_param_array[i][2][1]*11.60452501)
        Error_T_array.append(main_param_array[i][2][2]*11.60452501)
        Ne_array.append(main_param_array[i][3][1])
        Error_Ne_array.append(main_param_array[i][3][2])
        Mg_array.append(main_param_array[i][4][1])
        Error_Mg_array.append(main_param_array[i][4][2])
        Si_array.append(main_param_array[i][5][1])
        Error_Si_array.append(main_param_array[i][5][2])
        S_array.append(main_param_array[i][6][1])
        Error_S_array.append(main_param_array[i][6][2])
        Ca_array.append(main_param_array[i][7][1])
        Error_Ca_array.append(main_param_array[i][7][2])
        Fe_array.append(main_param_array[i][8][1])
        Error_Fe_array.append(main_param_array[i][8][2])
        norm_array.append(main_param_array[i][9][1])
        Error_norm_array.append(main_param_array[i][9][2])
        red_chi_sq_array.append(main_param_array[i][10][1])
        daxss_x123_slow_count.append(daxsslevel1['X123_SLOW_CORRECTED'][main_param_array[i][0]])

    x123_slow_count_avg = (daxsslevel1['X123_SLOW_CORRECTED'][start_index] + daxsslevel1['X123_SLOW_CORRECTED'][start_index+1] + daxsslevel1['X123_SLOW_CORRECTED'][start_index+2] + daxsslevel1['X123_SLOW_CORRECTED'][start_index+4])/4
    if(subAverage == 1):
        daxss_x123_slow_count_plot = daxss_x123_slow_count - x123_slow_count_avg
        for i in range (0,len(daxss_x123_slow_count_plot),1):
            if(daxss_x123_slow_count_plot[i]<0):
                daxss_x123_slow_count_plot[i]=0
    elif(subAverage == 0):
        daxss_x123_slow_count_plot = daxss_x123_slow_count
    daxss_time = date2num(time_str_array)

    fig = plt.figure(num=None, figsize=(12.4, 9.6), facecolor='w', edgecolor='k')
    ax1 = fig.add_axes([0.1, 0.4, 0.8, 0.55])
    plt.title("Temperature Fit",fontsize=12)
    ax1.xaxis.set_visible(False)
    ax1.plot_date(daxss_time, daxss_x123_slow_count_plot , 'o--', color='red', label="DAXSS Slow Counts")
    ax1.set_ylabel("DAXSS Slow CPS (Counts/sec)", color="red", fontsize=12)
    ax2 = ax1.twinx()
    ax2.errorbar(daxss_time, T_array, yerr=Error_T_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8, color='blue', label="Temp (MK)")
    ax2.set_ylabel("Temperature (MK)", color="blue", fontsize=12)
    if(subAverage == 1):
        plt.axhline(avg_fit_param_array[0][1]*11.60452501, linestyle='dashed', color='blue')
    plt.xlim([daxss_time[0], daxss_time[-1]])
    ax1.xaxis.set_major_formatter(DateFormatter('%m-%d-%H-%M-%S'))
    ax1.tick_params(rotation=45)
    plt.legend()

    ax3 = fig.add_axes([0.1, 0.15, 0.8, 0.25])
    plt.axhline(1, linestyle='dashed', color='grey')
    ax3.plot_date(daxss_time, red_chi_sq_array, 'o--', color='black', label="Reduced Chi Squared")
    plt.xlim([daxss_time[0], daxss_time[-1]])
    ax3.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M:%S'))
    ax3.tick_params(rotation=45)
    plt.xlabel('Date Time (MM-DD hh:mm:ss) UTC',fontsize=12)
    plt.ylabel('Reduced Chi Square', fontsize=12)
    plt.savefig(plots_path+'/' + "Flare_Temp_FIT.png")

    fig = plt.figure(num=None, figsize=(12.4, 9.6), facecolor='w', edgecolor='k')
    ax1 = fig.add_axes([0.1, 0.4, 0.8, 0.55])
    plt.title("Abundance Factor Fit",fontsize=12)
    ax1.xaxis.set_visible(False)
    ax1.plot_date(daxss_time, daxss_x123_slow_count_plot, 'o--', color='red', label="DAXSS Slow Counts")
    ax1.set_ylabel("DAXSS Slow CPS (Counts/sec)", color="red", fontsize=12)
    ax2 = ax1.twinx()
    #ax2.errorbar(daxss_time, Ne_array, yerr=Error_Ne_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8, color='blue', label="Ne AF")
    ax2.errorbar(daxss_time, Mg_array, yerr=Error_Mg_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8,  color='green', label="Mg AF")
    ax2.errorbar(daxss_time, Si_array, yerr=Error_Si_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8, color='black', label="Si AF")
    ax2.errorbar(daxss_time, S_array,  yerr=Error_S_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8, color='cyan', label="S AF")
    ax2.errorbar(daxss_time, Ca_array, yerr=Error_Ca_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8, color='yellow', label="Ca AF")
    ax2.errorbar(daxss_time, Fe_array, yerr=Error_Fe_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8, color='magenta', label="Fe AF")

    if (subAverage == 1):
        #plt.axhline(avg_fit_param_array[1][1], linestyle='dashed', color='blue')
        plt.axhline(avg_fit_param_array[2][1], linestyle='dashed', color='green')
        plt.axhline(avg_fit_param_array[3][1], linestyle='dashed', color='black')
        plt.axhline(avg_fit_param_array[4][1], linestyle='dashed', color='cyan')
        plt.axhline(avg_fit_param_array[5][1], linestyle='dashed', color='yellow')
        plt.axhline(avg_fit_param_array[6][1], linestyle='dashed', color='magenta')

    ax2.set_ylabel("Abundance Factor", color="black", fontsize=12)
    plt.xlim([daxss_time[0], daxss_time[-1]])
    ax1.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M:%S'))
    ax1.tick_params(rotation=45)
    plt.legend()
    ax3 = fig.add_axes([0.1, 0.15, 0.8, 0.25])
    plt.axhline(1, linestyle='dashed', color='black')
    ax3.plot_date(daxss_time, red_chi_sq_array, 'o--', color='black', label="Reduced Chi Squared")
    plt.xlim([daxss_time[0], daxss_time[-1]])
    ax3.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M:%S'))
    ax3.tick_params(rotation=45)
    plt.xlabel('Date Time (MM-DD hh:mm:ss) UTC', fontsize=12)
    plt.ylabel('Reduced Chi Square',color='black', fontsize=12)
    plt.savefig(plots_path +'/' + "Flare_AF_FIT.png")

    fig = plt.figure(num=None, figsize=(12.4, 9.6), facecolor='w', edgecolor='k')
    ax1 = fig.add_axes([0.1, 0.4, 0.8, 0.55])
    plt.title("Norm Fit", fontsize=12)
    ax1.xaxis.set_visible(False)
    ax1.plot_date(daxss_time, daxss_x123_slow_count_plot, 'o--', color='red', label="DAXSS Slow Counts")
    ax1.set_ylabel("DAXSS Slow CPS (Counts/sec)", color="red", fontsize=12)
    ax2 = ax1.twinx()
    ax2.errorbar(daxss_time, norm_array, yerr=Error_norm_array, fmt='o--', ms=0.5, capsize=1.0, lw=0.8, color='blue', label="Norm")
    ax2.set_ylabel("Norm", color="blue", fontsize=12)
    if (subAverage == 1):
        plt.axhline(avg_fit_param_array[7][1], linestyle='dashed', color='blue')
    plt.xlim([daxss_time[0], daxss_time[-1]])
    ax1.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M:%S'))
    ax1.tick_params(rotation=45)
    plt.legend()

    ax3 = fig.add_axes([0.1, 0.15, 0.8, 0.25])
    plt.axhline(1, linestyle='dashed', color='black')
    ax3.plot_date(daxss_time, red_chi_sq_array, 'o--', color='black', label="Reduced Chi Squared")
    plt.xlim([daxss_time[0], daxss_time[-1]])
    ax3.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M:%S'))
    ax3.tick_params(rotation=45)
    plt.xlabel('Date Time (MM-DD hh:mm:ss) UTC', fontsize=12)
    plt.ylabel('Reduced Chi Square',color='black', fontsize=12)
    plt.savefig(plots_path+'/' + "Flare_Norm_FIT.png")
    plt.show()


# Primary HDU - Header
hdr_dummy = fits.Header()
hdr_data = fits.Header()
hdr_dummy['MISSION'] = "InspireSat-1"
hdr_dummy['TELESCOP'] = "InspireSat-1"
hdr_dummy['INSTRUME'] = "DAXSS"
hdr_dummy['ORIGIN'] = "LASP"
hdr_dummy['CREATOR'] = "DAXSSPlotterUtility_v1"
hdr_dummy['CONTENT'] = "Type-I PHA file"

# Data Header
hdr_data['MISSION'] = "InspireSat-1"
hdr_data['TELESCOP'] = "InspireSat-1"
hdr_data['INSTRUME'] = "DAXSS"
hdr_data['ORIGIN'] = "LASP"
hdr_data['CREATOR'] = "DAXSSPlotterUtility_v1"
hdr_data['CONTENT'] = "SPECTRUM"
hdr_data['HDUCLASS'] = "OGIP"
hdr_data['LONGSTRN'] = "OGIP 1.0"
hdr_data['HDUCLAS1'] = "SPECTRUM"
hdr_data['HDUVERS1'] = "1.2.1"
hdr_data['HDUVERS'] = "1.2.1"

hdr_data['AREASCAL'] = "1"
hdr_data['BACKSCAL'] = "1"
hdr_data['CORRSCAL'] = "1"
hdr_data['BACKFILE'] = "none"

hdr_data['RESPFILE'] = "FITS_Files/minxss_fm3_RMF.fits"
hdr_data['ANCRFILE'] = "FITS_Files/minxss_fm3_ARF.fits"

hdr_data['CHANTYPE'] = "PHA"
hdr_data['POISSERR'] = "F"

hdr_data['CORRFILE'] = "none"
hdr_data['EXTNAME'] = 'SPECTRUM'
hdr_data['FILTER'] = "Be/Kapton"
hdr_data['EXPOSURE'] = "9"
hdr_data['DETCHANS'] = "1000"
hdr_data['GROUPING'] = "0"


def main():
    global main_param_array, start_index, subtractAverage
    # Set Analysis Parameters
    doy_folder_name = 'DOY065'
    main_param_array = []
    start_index = 315
    #end_index = 4720
    end_index = 387
    subtractAverage = 1

    daxss_level1_file_path = "netCDF_Files/DAXSS/daxss_solarSXR_level1_2022-02-14-mission_v2.0.0.ncdf"
    daxsslevel1 = nc.Dataset(daxss_level1_file_path)

    if (subtractAverage == 1):
        folder_name = 'Without_Quiet'
    else:
        folder_name = 'With_Quiet'

    # Parent Directory path
    global plots_path, pha_path, logs_path
    parent_dir = "CaseStudyAnalysis"
    current_analysis_dir = parent_dir + '/' + doy_folder_name + '_' + folder_name
    plots_path = os.path.join(current_analysis_dir, 'Plot_Files')
    pha_path = os.path.join(current_analysis_dir, 'PHA_Files')
    logs_path = os.path.join(current_analysis_dir, 'LOG_Files')

    if os.path.exists(current_analysis_dir):
        shutil.rmtree(current_analysis_dir)
    os.makedirs(current_analysis_dir)

    os.mkdir(plots_path)
    os.mkdir(pha_path)
    os.mkdir(logs_path)

    if(subtractAverage == 1):
        avg_fit_param_array = generateAvgSpec(daxsslevel1, start_index)
        start_index += 5

    for i in range(start_index, end_index, 1):
        generateFITSFile(daxsslevel1, i, subtractAverage)
        print("Success Fit Number Complete", i)

    generateParamPlot(daxsslevel1, subtractAverage, avg_fit_param_array)
    print(main_param_array)

if __name__ == "__main__":
    main()

#main_param_array[i][9][1] > 5 or main_param_array[i][9][1] < 0.4