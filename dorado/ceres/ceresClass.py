import warnings
warnings.filterwarnings('ignore')
# import sys
# import os

# import numpy as np
import ccdproc
from astropy.time import Time
from astropy.table import QTable, Table
import astroalign as aa
from astropy.wcs import WCS
# from astropy.utils.console import ProgressBar, ProgressBarOrSpinner
from tqdm import tqdm
# from astropy.coordinates import SkyCoord as acoord
# import astropy.units as un
from astropy.io import fits

from astropy.nddata.ccddata import CCDData

# photometry imports
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.aperture import CircularAperture, aperture_photometry, CircularAnnulus
from photutils import DAOStarFinder
from astropy.stats import mad_std
# from astroquery.simbad import Simbad

from ..timeseries import timeSeries

'''
Ceres is the handler of image series in Dorado,
'''

__all__ = ['Ceres']

class Ceres:
    '''
        The Ceres class encapsulates a set of astronomical data from a single night of observation.
        Ceres can handle multiple stacks of data in different filters and perform a variety of
        actions on them in an orgainized fashion.

        Attributes
        ----------

        filters: dictionary

        data: Stack array

        bias: CCDdata

        time: 'astropy.Time' 

        datestr: str

        
    '''
    def __init__(self, filters = {}, data = [], bias = None, time = None, datestr = None):
        # metadata
        self.filters = filters
        self.data = data
        self.bias = bias
        self.time = time
        self.datestr = datestr
        
        # location, weather, timezone, camera, observer array

        # refering instance of clippy to call and save later?
        # or call clippy directly and feed it a ceres object

        try:
            self.date = Time(int(self.time.mjd), format = 'mjd')
        except:
            self.date = None
        
        if datestr == None:
            try:
                day = str(self.date.ymdhms['day'])
                day2 = str(self.date.ymdhms['day'] + 1)
                month = str(self.date.ymdhms['month'])

                if self.date.ymdhms['day'] < 10:
                    day = '0' + str(self.date.ymdhms['day'])
                    if self.date.ymdhms['day'] < 9:
                        day2 = '0' + str(self.date.ymdhms['day'] + 1)
                        
                if self.date.ymdhms['month'] < 10:
                    month = '0' + str(self.date.ymdhms['month'])

                self.datestr = str(self.date.ymdhms['year']) + '-' + month + '-' + day + '+' + day2
            except:
                self.datestr = datestr
        
    def add_stack(self, stack):
        # eventually stacks themelves should have some metadata 
        # to denote stuff like calibration status
        self.filters[stack.filter] = len(self.data)
        self.data.append(stack)

    def rem_stack(self, filter):
        del self.data[self.filters[filter]]
        # delete time strings

    def calibrate(self, filter):
        # for bla in series: add bias corrected = True to header
        stack = self.data[self.filters[filter]]
        flat = stack.flat
        bias = self.bias
        c_series = []
        # with ProgressBar(len(stack.data)) as bar:
        print('Calibrating')
        for im in tqdm(stack.data, colour = 'green'):
            # bar.update()
            im.data = im.data.astype('uint16') 
            flat.data = flat.data.astype('uint16') 
            bias.data = bias.data.astype('uint16') 
            im = ccdproc.ccd_process(im, master_bias = bias, master_flat = flat)
            im.data = im.data.astype('uint16') 
            c_series.append(im)
        self.data[self.filters[filter]].data = c_series
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

    def getWCS(self, filter, filer, alignto = None, cache = True):
        series = self.data[self.filters[filter]]
        if alignto == None:
            alignto = series.alignTo
        if cache:
            hdulist = fits.open(filer.dordir / 'cache' / 'astrometryNet' / 'solved.fits') 
            self.data[self.filters[filter]].wcs = WCS(hdulist[0].header, hdulist)
            self.data[self.filters[filter]].solved = CCDData.read(filer.dordir / 'cache' / 'astrometryNet' / 'solved.fits')
            hdulist.close()
        else:
            toalign = series.data[alignto]
            fname, cachedir = filer.mkcacheObj(toalign, 'astrometryNet')
            path = [cachedir, fname]
            writearray = [cachedir, 'solved.fits']
            solved, wcs_header = filer.plate_solve(path, writearray = writearray)
            filer.delcacheObj( fname, 'astrometryNet')
            self.data[self.filters[filter]].wcs = WCS(wcs_header)
            self.data[self.filters[filter]].solved = solved

    def align(self, filter, filer, alignto = None, getWCS = True, cache = False):
        series = self.data[self.filters[filter]]
        if alignto == None:
            alignto = series.alignTo
        toalign = series.data[alignto]
        ## TODO :: make this use ceres.getWCS()
        if getWCS:
            if cache:
                toalign =  CCDData.read(filer.dordir / 'cache' / 'astrometryNet' / 'solved.fits', unit = filer.unit)
                hdulist = fits.open(filer.dordir / 'cache' / 'astrometryNet' / 'solved.fits') 
                self.data[self.filters[filter]].wcs = WCS(hdulist[0].header, hdulist)
                hdulist.close()
                self.data[self.filters[filter]].solved = toalign
            else:
                fname, cachedir = filer.mkcacheObj(toalign, 'astrometryNet')
                path = [cachedir, fname]
                writearray = [cachedir, 'solved.fits']
                solved, wcs_header = filer.plate_solve(path, writearray = writearray)
                toalign = solved
                filer.delcacheObj( fname, 'astrometryNet')
                self.data[self.filters[filter]].wcs = WCS(wcs_header)
                self.data[self.filters[filter]].solved = solved
                # delete cache object
                # save solved to target

        aa_series = []
        skipped = []
        ## TODO :: fix this progressbar so it prints on one line then updates that line.
        # with ProgressBar(len(series.data)) as bar:
        print('Aligning')
        for image in tqdm(series.data, colour = 'green'):
            # bar.update()
            try:
                img, _ = aa.register(image.data, toalign.data)
                image.data = img
                aa_series.append(image)
            except:
                skipped.append(image)
                # print('Image skipped')
        if len(skipped) != 0:
            print(len(skipped), ' images skipped.')
        self.data[self.filters[filter]].data = aa_series
        self.data[self.filters[filter]].aligned = True

    def dorphot(self, filter, toi, control_toi = None, shape = 21, unc = 0.1):
        # get seeing from PSF
        stack = self.data[self.filters[filter]]
        # if no wcs, complain alot
        w = stack.wcs

        xy = w.wcs_world2pix(toi.coords.ra.deg, toi.coords.dec.deg, 1)
        ra = toi.coords.ra.deg
        dec = toi.coords.dec.deg
        # pos = Table(names=['x_0', 'y_0'], data = ([float(xy[0])], [float(xy[1])]))
        pos = [(float(xy[0]), float(xy[1]))]
        aperture = CircularAperture(pos, r = shape)
        annulus_aperture = CircularAnnulus(pos, r_in = shape + 2, r_out = shape + 5)
        apers = [aperture, annulus_aperture]

        if control_toi != None:
            xyc = w.wcs_world2pix(control_toi.coords.ra.deg, control_toi.coords.dec.deg, 1)
            # posc = Table(names=['x_0', 'y_0'], data = ([float(xyc[0])], [float(xyc[1])]))
            posc = [(float(xyc[0]), float(xyc[1]))]
            aperturec = CircularAperture(posc, r = shape)
            annulus_aperturec = CircularAnnulus(posc, r_in = shape + 2, r_out = shape + 5)
            apersc = [aperturec, annulus_aperturec]


        times = []
        exptimes = []
        ray = []
        decx = []
        x = []
        y = []
        flux = []
        fluxunc = []
        apsum = []
        apsum_unc = []


        print('Performing photometry')
        for image in tqdm(stack.data, colour = 'green'):
            error = unc * image.data
            results = aperture_photometry(image, apers, error = error)
            bkg_mean = results['aperture_sum_1'] / annulus_aperture.area
            bkg_sum = bkg_mean * aperture.area
            results['flux_fit'] = results['aperture_sum_0'] - bkg_sum
            
            times.append(Time(image.header['DATE-OBS']))
            exptimes.append(image.header['EXPTIME'])
            ray.append(ra)
            decx.append(dec)
            x.append(results['xcenter'][0])
            y.append(results['ycenter'][0])
            # x.append(results['x_fit'][0])
            # y.append(results['y_fit'][0])

            if control_toi != None:
                resultsc = aperture_photometry(image, apersc, error = error)
                bkg_meanc = resultsc['aperture_sum_1'] / annulus_aperturec.area
                bkg_sumc = bkg_meanc * aperturec.area
                resultsc['flux_fit'] = resultsc['aperture_sum_0'] - bkg_sumc

                apsum.append(results['flux_fit'][0] - resultsc['flux_fit'][0])
                flux.append((results['flux_fit'][0] - resultsc['flux_fit'][0])/image.header['EXPTIME'])
            else:
                apsum.append(results['flux_fit'][0])
                flux.append(results['flux_fit'][0]/image.header['EXPTIME'])

            fluxunc.append(results['aperture_sum_err'][0]) ## TODO:: modify this to account for exposure time and control
            apsum_unc.append(results['aperture_sum_err'][0])

        ts = timeSeries(times = times, flux = flux, exptimes = exptimes, x = x, y = y, ra = ray, dec = decx, flux_unc = fluxunc, apsum = apsum, apsum_unc = apsum_unc)
        toi.filters[filter] = len(toi.ts)
        toi.ts.append(ts)

