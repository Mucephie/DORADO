# exotest.py
import dorado
import astroalign as aa
import datetime
import numpy as np
from astroquery.simbad import Simbad
from astropy.nddata import CCDData
from astropy.table import QTable, Table
from astropy.io import fits
from astropy import units as un

# time imports
from astropy.time import Time
from astropy.timeseries import TimeSeries

# photometry imports
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils import DAOStarFinder
from astropy.stats import mad_std

# coordinate imports
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord as acoord
from astropy.coordinates import Distance


## One might be interested in how long this all takes
START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)
print('\nFinding data...')





## Organize and prepare your data
directory = 'C:/Data/2020-22-04+05/'
# You can name your target object
target = 'WASP33b'
s = Simbad()
r = s.query_object(target)
r.pprint()
print(r.colnames)
co = acoord(ra = r['RA'], dec = r['DEC'], unit = (un.hourangle, un.deg), frame = 'icrs')

# Create a working directory
dorado.mkwrkdir(directory)


# Catalogue input data from the data directory
_, _, _, lights_list = dorado.checkdir(directory + 'lights/')
_, flats_list, _, _= dorado.checkdir(directory + 'flats/')
bias_list, _, _, _ = dorado.checkdir(directory + 'bias/')



## Calibrate your data
# Read the data in
# data_series = dorado.get_series(directory + 'lights/', lights_list)
# flats_series = dorado.get_series(directory, flats_list)
# bias_series = dorado.get_series(directory, bias_list)

# Produce master reduction images
bias = dorado.mastBias(directory + 'bias/', bias_list, caldir = directory)
flat = dorado.mastFlat(directory + 'flats/', flats_list, bias, caldir = directory)

# Calibrate data series, files will be written automatically
series = dorado.reduce_series(directory + 'lights/', lights_list, flat, bias, target, caldir = directory)

CAL_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', CAL_END_DATE_TIME)
print("Time elapsed for calibration: ", (CAL_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')





## Plate solve your data series
# Pick image to solve and store its data
psimg = series[0]
# Set up a file string for dorado to locate image
expath = '/wrk/calibrated/' 
im = target +  '-0_c.fit'
# image_file_path = dorado.file_string(directory + expath, im)
image_file_path = dorado.file_string(directory, im)

# Pass the data and file string to dorado
solved = dorado.plate_solve(psimg, image_file_path, write_fits = True, write_name = target + '_solved')
w = WCS(solved.header)
PLATE_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', PLATE_END_DATE_TIME)
print("Time elapsed for solving: ", (PLATE_END_DATE_TIME - CAL_END_DATE_TIME))
print("Time elapsed so far: ", (PLATE_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')





# Align each image with the plate solved image
aa_series = []
for image in series:
    img_r, footprint_r = aa.register(image.data, solved.data)
    image.data = img_r
    aa_series.append(image)

dorado.write_series(directory + '/wrk/aligned/', aa_series, target, '_a')
# You now have a fully calibrated and aligned data set with WCS data
AL_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', AL_END_DATE_TIME)
print("Time elapsed for alignment: ", (AL_END_DATE_TIME - PLATE_END_DATE_TIME))
print("Time elapsed so far: ", (AL_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')





## Data analysis
# Run PSF photometry on image one
sigma_psf = 2.0
bkgrms = MADStdBackgroundRMS()
std = bkgrms(aa_series[0])
bkg_sigma = mad_std(aa_series[0]) 
daofind = DAOStarFinder(fwhm = 4., threshold = 3. * bkg_sigma)  
s = daofind(aa_series[0].data) 

daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma = sigma_psf)

photometry = IterativelySubtractedPSFPhotometry(finder = daofind, 
                        group_maker = daogroup, bkg_estimator = mmm_bkg,
                         psf_model = psf_model, fitter = LevMarLSQFitter(),
                          niters = 1, fitshape = (11, 11))
result_tab = photometry(image = aa_series[0])

sources = Table()
sources['ra'], sources['dec'] = w.wcs_pix2world(result_tab['x_fit'], result_tab['y_fit'], 1) 
sources['x'] = result_tab['x_fit']
sources['y'] = result_tab['y_fit']
sources['flux'] = result_tab['flux_fit']
sources['flux_unc'] = result_tab['flux_unc']

# clean sources
sources = sources[sources['flux'] > 2000]
sources['mag_inst'] = -2.5 * np.log10(sources['flux'] / 15)
sources['mag_inst_unc'] = sources['flux_unc'] / ((sources['flux'] / 15) * 2.30258509)
sources.write('sources_test.csv', overwrite = True)
# fix positions to run PSF photometry on remaining images
target_px_x, target_px_y = w.wcs_pix2world(co)
print('target px: ', target_px_x, target_px_y)
target_s = sources[np.abs(sources[np.abs(sources['x'] - target_px_x) < 4]['y'] - target_px_y) < 4]
print('target in sources: ', target_s)

psf_model.x_0.fixed = True
psf_model.y_0.fixed = True
pos = Table(names=['x_0', 'y_0'], data=[sources['x'], sources['y']])

# ts = QTable(names=('time', 'flux', 'flux_unc', 'dmag', 'dmag_std', 'dmag_unc'))
ts = QTable(names=('time', 'dmag', 'dmag_std', 'dmag_unc'))


S1_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', S1_END_DATE_TIME)
print("Time elapsed for photometry 1: ", (S1_END_DATE_TIME - AL_END_DATE_TIME))
print("Time elapsed so far: ", (S1_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')

for image in aa_series:
    time = Time(image.header['DATE-OBS'], format = 'fits', scale = 'utc')

    std = bkgrms(image)
    psf_model = IntegratedGaussianPRF(sigma = sigma_psf)
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True

    photometry = IterativelySubtractedPSFPhotometry(finder = daofind, 
                        group_maker = daogroup, bkg_estimator = mmm_bkg,
                         psf_model = psf_model, fitter = LevMarLSQFitter(),
                          niters = 1, fitshape = (11, 11))
    result_tab = photometry(image = image, init_guesses=pos)
    result_tab['mag'] = -2.5 * np.log10(result_tab['flux_fit'] / 15)
    result_tab['mag_unc'] = result_tab['flux_unc'] / ((result_tab['flux_fit'] / 15) * 2.30258509)

    target_idx = result_tab[np.abs(result_tab[np.abs(result_tab['x_fit'] - target_px_x) < 4]['y_fit'] - target_px_y) < 4]
    result_tab['dmag'] = target_idx['mag'] - result_tab['mag']
    dmag = np.mean(result_tab['dmag'])
    dmag_std  = np.std(result_tab['dmag'])
    dmag_unc = target_idx['mag_unc']
    ts.add_row((time, dmag, dmag_std, dmag_unc))

ts.show_in_browser(jsviewer = True)
ts.write('timeseries_test.csv', overwrite = True)


# Construct lightcurve

# Locate times of max light

# Call ephemeris to get O-C

## exit script, it took this long
END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', END_DATE_TIME)
print("Time elapsed for run: ", (END_DATE_TIME - START_DATE_TIME))




# Finished time:  2021-01-01 19:47:59.883863
# Time elapsed for alignment:  0:11:47.918621
# Time elapsed so far:  0:16:19.697864





# Traceback (most recent call last):
#   File "c:/Users/mucep/OneDrive/Documents/GitHub/DORADO/dorado/examples/exotest.py", line 153, in <module>
#     sources['mag_inst'] = -2.5 * np.log10(sources['flux']/ aa_series[0].header['EXPTIME'])
#   File "C:\Users\mucep\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.8_qbz5n2kfra8p0\LocalCache\local-packages\Python38\site-packages\astropy\io\fits\header.py", line 148, in __getitem__
#     card = self._cards[self._cardindex(key)]
#   File "C:\Users\mucep\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.8_qbz5n2kfra8p0\LocalCache\local-packages\Python38\site-packages\astropy\io\fits\header.py", line 1708, in _cardindex
#     raise KeyError(f"Keyword {keyword!r} not found.")
# KeyError: "Keyword 'EXPTIME' not found."