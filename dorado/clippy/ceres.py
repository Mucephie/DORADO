import warnings
warnings.filterwarnings('ignore')

import numpy as np
import ccdproc
from astropy.time import Time
from astropy.timeseries import BinnedTimeSeries
from astropy.table import QTable, Table
import astroalign as aa
from astropy.wcs import WCS
from astropy.utils.console import ProgressBar
from astropy.coordinates import SkyCoord as acoord
import astropy.units as un

# photometry imports
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils import DAOStarFinder
from astropy.stats import mad_std
from astroquery.simbad import Simbad

'''
Ceres is the handler of series in Dorado,
'''

# __all__ = [Ceres]

class Ceres:

    def __init__(self):
        # calibrate
        # metadata
        self.filters = {}
        self.data = []
        self.bias = None
        self.flats = {}
        self.darks = {}
        self.time = {}
        # self.flats = {}
        # self.darks = {}
        self.time = None
        
        # date, location, weather, timezone, camera, observer array
        # refering instance of clippy to call and save later?
        # or call clippy directly and feed it a ceres object
        

        try:
            self.date = Time(int(self.time.mjd), format = 'mjd')
        except:
            self.date = None
        
    def add_stack(self, stack):
        # eventually stacks themelves should have some metadata 
        # to denote stuff like calibration status
        self.filters[stack.filter] = len(self.data)
        self.data.append(stack)
        # this should also extract the time strings

    def rem_stack(self, filter):
        del self.data[self.filters[filter]]
        # delete time strings

    def calibrate(self, filter):
        # for bla in series: add bias corrected = True to header
        flat = self.flats[filter]
        proc = self.data[self.filters[filter]]
        flat = proc.flat
        bias = self.bias
        proc = self.data[filter]
        for p in range(len(proc)):
            proc[p] = ccdproc.ccd_process(proc[p], master_bias = bias, master_flat = flat)
        self.data[filter] = proc
    
        with ProgressBar(len(proc.data)) as bar:
            for p in range(len(proc.data)):
                bar.update()
                proc.data[p] = ccdproc.ccd_process(proc.data[p], master_bias = bias, master_flat = flat)
        self.data[self.filters[filter]].data = proc
        self.data[self.filters[filter]].calibrated = True

    def imarith(self, filter, operator, operand):
        # mod to check datatype using type()
        # mod to remove im_count and make possible to use single image
        # mod to accomodate CCDdata object
        series = self.data[self.filters[filter]]
        for i in range(len(series)):
            if (operator == '+'):
                series[i].data = series[i].data  + operand
            elif (operator == '-'):
                series[i].data = series[i].data - operand
            elif (operator == '/'):
                series[i].data = series[i].data  / operand
            elif (operator == '*'):
                series[i].data = series[i].data  * operand
        
        self.data[self.filters[filter]] = series

    def align(self, filter, clippy, alignto = 0, getWCS = True):
        series = self.data[self.filters[filter]]
        toalign = series.data[alignto]
        if getWCS:
            fname, cachedir = clippy.mkcacheObj(toalign, 'astrometryNet')
            path = [cachedir, fname]
            writearray = [cachedir, 'solved.fits']
            solved, wcs_header = clippy.plate_solve(path, writearray = writearray)
            toalign = solved
            self.data[self.filters[filter]].wcs = WCS(wcs_header)
            # delete cache object
            # save solved to target

        aa_series = []
        with ProgressBar(len(series.data)) as bar:
            for image in series.data:
                bar.update()
                img, _ = aa.register(image.data, toalign.data)
                aa_series.append(img)
        self.data[self.filters[filter]].data = aa_series
        self.data[self.filters[filter]].aligned = True

    # save to wrk

    def dorphot(self, filter, zellars):
        # get seeing from PSF
        stack = self.data[self.filters[filter]]
        # if no wcs, complain alot
        w = stack.wcs
        print(w.wcs_world2pix(zellars.coords.ra.deg, zellars.coords.dec.deg, 1))
        xy = w.wcs_world2pix(zellars.coords.ra.deg, zellars.coords.dec.deg, 1)
        pos = Table(names=['x_0', 'y_0'], data=[float(xy[0]), float(xy[1])])
        sigma_psf = 2.0
        bkg_sigma = mad_std(stack.data[0]) 
        daofind = DAOStarFinder(fwhm = 4., threshold = 3. * bkg_sigma)  
        daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
        mmm_bkg = MMMBackground()

        psf_model = IntegratedGaussianPRF(sigma = sigma_psf)
        psf_model.x_0.fixed = True
        psf_model.y_0.fixed = True

        fitter = LevMarLSQFitter()
        bkgrms = MADStdBackgroundRMS()

        ts = Table([[], [], [], [], [], [], [], []], names=('x', 'y', 'ra', 'dec', 'flux', 'flux_unc', 'time', 'exptime'), meta={'name': filter})

        # if radec get xy
        with ProgressBar(len(stack.data)) as bar:
            for image in stack.data:
                bar.update()
                photometry = IterativelySubtractedPSFPhotometry(finder = daofind, group_maker = daogroup, bkg_estimator = mmm_bkg,
                        psf_model = psf_model, fitter = LevMarLSQFitter(), niters = 1, fitshape = (21, 21))

                results = photometry(image = image, init_guesses= pos)


                time = Time(image.header['DATE-OBS'], format='fits')
                exptime = image.header['EXPTIME']
                [ra, dec] = w.wcs_pix2world(results['x_fit'], results['y_fit'])

                ts.add_row(results['x_fit'], results['y_fit'], ra, dec, results['flux_fit'], results['flux_unc'], time, exptime)
                # instrumental magnitude
        zellars.filters[filter] = len(zellars.ts)
        zellars.ts.append(ts)



class zellars:
    def __init__(self, name):
        # get it because zellars is the canadian target?
        self.name = name
        s = Simbad()
        r = s.query_object(self.name)
        # r.pprint()
        # print(r.colnames)
        self.coords = acoord(ra = r['RA'], dec = r['DEC'], unit = (un.hourangle, un.deg), frame = 'icrs')
        


class Stack:
    def __init__(self, data, filter = '', times = [], calibrated = None, aligned = None, target = ''):
        self.data = data
        self.filter = filter
        self.length = len(data)
        self.calibrated = calibrated
        self.aligned = aligned
        self.target = target
        self.target_info = {} # dictionary of values
        self.times = times # put together a check for if filled, if not, try to find it.

        # include things like flux uncertainty etc.
        # include wcs


        if filter == '':
            try:
                filter = data[0].header['filter']
            except:
                filter = ''
        