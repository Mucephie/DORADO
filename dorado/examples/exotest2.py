from astropy.utils.misc import silence
silence()
import warnings
warnings.filterwarnings('ignore')

# exotest.py
import dorado
import astroalign as aa
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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


# read in test data
hdulist = fits.open(directory + 'new-image.fits') # 'wcs.fits')
# print(hdulist)
w = WCS(hdulist[0].header, hdulist)
hdulist.close()
# print(co.ra, co.dec)
print('coord test: ', w.wcs_world2pix(co.ra.deg, co.dec.deg, 1) )


# Catalogue input data from the data directory
_, _, _, lights_list = dorado.checkdir(directory + '/wrk/aligned/')


aa_series = dorado.get_series(directory + '/wrk/aligned/', lights_list)


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
target_px_x, target_px_y = w.wcs_world2pix(co.ra.deg, co.dec.deg, 1)
print('\n\n\ntarget px: ', target_px_x, target_px_y)
target_s_temp = sources[np.abs(sources['x'] - target_px_x) < 4]
target_s = target_s_temp[np.abs(target_s_temp['y'] - target_px_y) < 4]
print('target in sources: ')
target_s.pprint()

psf_model.x_0.fixed = True
psf_model.y_0.fixed = True
pos = Table(names=['x_0', 'y_0'], data=[sources['x'], sources['y']])

# ts = QTable(names=('time', 'flux', 'flux_unc', 'dmag', 'dmag_std', 'dmag_unc'))
ts = QTable(names=('time', 'dmag', 'dmag_std', 'dmag_unc'))


S1_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', S1_END_DATE_TIME)
print('Finished loading and photometry 1')
print("Time elapsed so far: ", (S1_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')

for idx in range(len(aa_series[1::])):
    image = aa_series[idx + 1]
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
    # result_tab = result_tab[result_tab['flux_fit'] > 2000]
    result_tab['mag'] = -2.5 * np.log10(np.array(result_tab['flux_fit']) / 15)
    result_tab['mag_unc'] = np.array(result_tab['flux_unc']) / ((np.array(result_tab['flux_fit']) / 15) * 2.30258509)

    # target_idx = result_tab[np.abs(result_tab[np.abs(result_tab['x_fit'] - target_px_x) < 4]['y_fit'] - target_px_y) < 4]
    target_idx_temp = result_tab[np.abs(result_tab['x_fit'] - target_px_x) < 4]
    target_idx = target_idx_temp[np.abs(target_idx_temp['y_fit'] - target_px_y) < 4]
    
    result_tab['dmag'] = target_idx['mag'] - result_tab['mag']
    dmag = np.mean(result_tab['dmag'])
    dmag_std  = np.std(result_tab['dmag'])
    dmag_unc = target_idx['mag_unc']
    print('row ', idx, ': ', time.fits, float(dmag), float(dmag_std), float(dmag_unc))
    ts.add_row((float(time.mjd), float(dmag), float(dmag_std), float(dmag_unc)))

ts.show_in_browser(jsviewer = True)
ts.write('timeseries_test.csv', overwrite = True)


# Construct lightcurve
plt.figure()

plt.scatter(ts['time'], ts['dmag']) 

plt.show()
# Locate times of max light

# Call ephemeris to get O-C

## exit script, it took this long
END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', END_DATE_TIME)
print("Time elapsed for run: ", (END_DATE_TIME - START_DATE_TIME))
