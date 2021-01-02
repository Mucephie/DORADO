from astroquery.simbad import Simbad
import numpy as np
from astropy import units as un
from astropy.coordinates import SkyCoord as acoord
from astropy.coordinates import Distance
from astropy.table import QTable, Table
from photutils.detection import IRAFStarFinder
from photutils import DAOStarFinder
from astropy.stats import mad_std
from astropy.wcs import WCS

from astropy.nddata import CCDData
from astropy.io import fits
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IterativelySubtractedPSFPhotometry

directory = 'E:/DORADO/2020-22-04+05/'

target = 'WASP33b'
s = Simbad()
r = s.query_object(target)
r.pprint()
print(r.colnames)
# print(r['RA'], r['DEC'])


class toi:
    # work on a class for targets of interest

    # initialization of class
    def __init__(self, name):
        # identification
        self.name = name # toi name
        self.id = None # toi ID
        self.object_type = None # object type of toi
        # celestial coordinates
        self.ra = None # right accension of toi
        self.dec = None # declination of toi
        self.distance = None # distance quantity of toi
        self.z = None # redshift of toi
        self.co = None # coordinate object of toi
        # initialize photometry
        self.flux = QTable(names=('filter', 'flux', 'unc'))
        self.mag = QTable(names=('filter', 'mag', 'unc'))
        # image coordinates 
        self.xpx = None # x coordinate of toi
        self.ypx = None # y coordinate of toi

    # class functions
    def add_flux(self, filt, flux, unc, add_mag = False):
        self.flux.add_row((filt, flux, unc))
        if add_mag:
            mag_unc = unc / (flux * 2.30258509)
            mag = -2.5 * np.log10(flux)
            self.add_mag(filt, mag, mag_unc)
                
    def add_mag(self, filt, mag, unc):
        self.mag.add_row((filt, mag, unc))

    def init_co(self, ra, dec, distance = None, z = None, frame = 'icrs', units = None):
        # initialize toi coordinates
        if distance != None:
            self.distance = Distance(distance)
            if units != None:
                self.co = acoord(ra = ra, dec = dec, frame = frame, distance = distance, unit = units)
            else:
                self.co = acoord(ra = ra, dec = dec, frame = frame, distance = distance)

            self.ra = self.co.ra
            self.dec = self.co.dec
            if z != None:
                self.z = z
            else:
                self.z = self.distance.z

        elif z != None:
            self.distance = Distance(z = z)
            if units != None:
                self.co = acoord(ra = ra, dec = dec, frame = frame, distance = distance, unit = units)
            else:
                self.co = acoord(ra = ra, dec = dec, frame = frame, distance = distance) 
            
            self.ra = self.co.ra
            self.dec = self.co.dec
            self.z = z

        else: 
            if units != None:
                self.co = acoord(ra = ra, dec = dec, frame = frame, unit = units)
            else:
                self.co = acoord(ra = ra, dec = dec, frame = frame)
            self.ra = self.co.ra
            self.dec = self.co.dec




wasp = toi(target)
wasp.init_co(ra = r['RA'], dec = r['DEC'], units = (un.hourangle, un.deg))
wasp.object_type = 'STAR'
# print(wasp.co.to_string('hmsdms'))


# read in test data
hdulist = fits.open(directory + 'wcs.fits')
# print(hdulist)
w = WCS(hdulist[0].header, hdulist)
hdulist.close()
axy = Table.read(directory + 'axy.fits')
image  = CCDData.read(directory + 'new-image.fits', unit = un.adu)

# set up PSF stuff
sigma_psf = 2.0
bkgrms = MADStdBackgroundRMS()
std = bkgrms(image)
iraffind = IRAFStarFinder(threshold = 3.5 * std, fwhm = sigma_psf * gaussian_sigma_to_fwhm,
                              minsep_fwhm = 0.01, roundhi = 5.0, roundlo = -5.0,
                               sharplo = 0.0, sharphi = 2.0)

bkg_sigma = mad_std(image) 
daofind = DAOStarFinder(fwhm = 4., threshold = 3. * bkg_sigma)  
s = daofind(image.data) 
# s.show_in_browser(jsviewer = True)

daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma = sigma_psf)

# perform photometry
photometry = IterativelySubtractedPSFPhotometry(finder = daofind, 
                        group_maker = daogroup, bkg_estimator = mmm_bkg,
                         psf_model = psf_model, fitter = LevMarLSQFitter(),
                          niters = 1, fitshape = (11, 11))
# photometry = IterativelySubtractedPSFPhotometry(finder = iraffind, 
#                         group_maker = daogroup, bkg_estimator = mmm_bkg,
#                          psf_model = psf_model, fitter = LevMarLSQFitter(),
#                           niters = 1, fitshape = (11, 11))
result_tab = photometry(image = image)

sources = Table()
sources['ra'], sources['dec'] = w.wcs_pix2world(result_tab['x_fit'], result_tab['y_fit'], 1) 
sources['x'] = result_tab['x_fit']
sources['y'] = result_tab['y_fit']
sources['flux'] = result_tab['flux_fit']
sources['flux_unc'] = result_tab['flux_unc']

# clean sources
sources = sources[sources['flux'] > 2000]
sources['mag_inst'] = -2.5 * np.log10(result_tab['flux_fit'])
sources['mag_inst_unc'] = sources['flux_unc'] / (sources['flux'] * 2.30258509)
sources.show_in_browser(jsviewer = True)


# wasp33 B-V 0.297
# V mag 8.14

# dorado.plt_stars(image, sources['x'], sources['x'], rs)