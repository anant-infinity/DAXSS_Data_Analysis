from astropy.io import fits
import numpy as np
import netCDF4 as nc
#fits_image_filename = daxss_level1_file_path = 'C:/Users/278an/Desktop/science_analysis/301_Data/PYXSPEC/minxss_fm3_RMF_ARF.fits'
#hdul = fits.open(fits_image_filename)
#print(hdul.info())

#n = np.arange(100.0)

# Path to Level-1 File
daxss_level1_file_path = 'C:/Users/278an/Desktop/science_analysis/301_Data/DAXSS/level1/daxss_solarSXR_level1_2022-02-14-mission_v1.0.0_updated.ncdf'
# Import File as a netCDF Dataset
daxsslevel1 = nc.Dataset(daxss_level1_file_path)
index_quiet_sun = 0


# Creating the FITS file for DAXSS
# Primary HDU - Header
hdr = fits.Header()
hdr['EXTNAME'] = 'SPECTRUM'
hdr['TELESCOP'] = "InspireSat-1"
hdr['INSTRUME'] = "DAXSS"
hdr['INSTRUME'] = "DAXSS"
hdr['FILTER'] = "Be/Kapton"
hdr['EXPOSURE'] = "9 seconds"
hdr['BACKFILE'] = "NA"
hdr['BACKSCAL'] = "NA"
hdr['CORRFILE'] = "NA"
hdr['RESPFILE'] = "minxss_fm3_RMF_ARF.fits"
hdr['ANCRFILE'] = "minxss_fm3_RMF_ARF.fits"
hdr['AREASCAL'] = "NA"
hdr['HDUCLASS'] = "OGIP"
hdr['HDUCLAS1'] = "SPECTRUM"
hdr['HDUVERS'] = "1.2.1"
hdr['POISSERR'] = "NA"
hdr['CHANTYPE'] = "PHA"
hdr['DETCHANS'] = "1000"

dummy_primary = fits.PrimaryHDU()

c1 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]
c2 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]
c3 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]
c4 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]
c5 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]
c6 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]
c7 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]
c8 = daxsslevel1['IRRADIANCE'][index_quiet_sun,:]

# Data
hdu_data = fits.BinTableHDU.from_columns(
    [fits.Column(name='CHANNEL', format='20A', array=c1),
     fits.Column(name='RATE', format='20A', array=c2),
     fits.Column(name='STAT_ERR', format='20A', array=c3),
     fits.Column(name='SYS_ERR', format='20A', array=c4),
     fits.Column(name='QUALITY', format='20A', array=c5),
     fits.Column(name='GROUPING', format='20A', array=c6),
     fits.Column(name='AREASCAL', format='20A', array=c7),
     fits.Column(name='BACKSCAL', format='E', array=c8)],header=hdr)
# Creating and Storing the FITS File
hdul = fits.HDUList([dummy_primary,hdu_data])
hdul.writeto('DAXSS_Spectrum2.fits',overwrite=True)