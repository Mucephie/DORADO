from tqdm import tqdm
import numpy as np

from ..core.coreClass import *
from ..timeseries.timeseriesClass import timeSeries
import ccdprocx

from astropy.time import Time
from astropy.nddata.ccddata import CCDData
import astroalign as aa
from astropy.io import fits
from astropy.wcs import WCS

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

from photutils.aperture import CircularAperture, aperture_photometry, CircularAnnulus


__all__ = ['aicoPhot']

class aicoPhot:

    def __init__(self):
        # TODO make observatory class
        self.temp = None
        # list of calibration frames from disk
        

    
    def calibrate(self, cr,  filter, use_med_cr_removal = False, rb = 0):
        '''
        calibrate performs 'pre-processing' CCDData calibration to a data stack within a ceres object. This 
        involves Bias and Flatfield correction and optional median value based cosmic ray removal based on the
        ccdprocx.cosmicray_median() method.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to calibrate.
        filter: str
            String representation of the relevent filter.
        use_med_cr_removal: boolean
            controls whether to use cosmic ray removal in calibration, may affect runtime. Default is False.
        rb: int
            the rbox value for ccdprocx.cosmicray_median(). Default is 0.
        
        '''
        # for bla in series: add bias corrected = True to header
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        flat = stack.flat
        bias = Dorado.ceres[Dorado.ceres_keys[cr]].bias
        c_series = []

        print('Calibrating')
        for im in tqdm(stack.data, colour = 'green'):
            # bar.update()
            im.data = im.data.astype('uint16') 
            flat.data = flat.data.astype('uint16') 
            bias.data = bias.data.astype('uint16') 
            im = ccdprocx.ccd_process(im, master_bias = bias, master_flat = flat)
            if use_med_cr_removal:
                im = ccdprocx.cosmicray_median(im, rbox = rb)
            im.data = im.data.astype('uint16') 
            c_series.append(im)
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].data = c_series
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].calibrated = True
        
    def imarith(self, cr, filter, operator, operand):
        '''
        imarith is a basic replication of the IRAF image arithmatic tool imarith.
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform image arithmatic on.
        filter: str
            String representation of the relevent filter.
        operator: string
            String of operator to be used for arithmatic. Supported operations are '+', '-', '/', and '*'
        '''
        # TODO should this be in stack? like have a wrapper here?
        # mod to check datatype using type()
        # mod to remove im_count and make possible to use single image NOTE might already be done
        # mod to accomodate CCDdata object NOTE might already be done
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        for i in range(len(series)):
            if (operator == '+'):
                series[i].data = series[i].data  + operand
            elif (operator == '-'):
                series[i].data = series[i].data - operand
            elif (operator == '/'):
                series[i].data = series[i].data  / operand
            elif (operator == '*'):
                series[i].data = series[i].data  * operand
        
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]] = series
        

    def getWCS(self, cr, filter, alignto = None, cache = True):
        '''
        getWCS obtains WCS information for an image either via previously solved data in the cache or
        by passing the image to astrometryNet via dorado.plate_solve().
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to get WCS information for.
        filter: str
            String representation of the relevent filter.
        alignto: int
            Index of image to use as reference. Default is stack.alignto
        cache: boolean
            Controls whether to call astrometryNet for solve or use solved result stored in cache 
            from previous run. Default is True.
        '''
        # TODO mod so cache results are target specific
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        if alignto == None:
            alignto = series.alignTo
        if cache:
            hdulist = fits.open(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits') 
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(hdulist[0].header, hdulist)
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = CCDData.read(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits')
            hdulist.close()
        else:
            toalign = series.data[alignto]
            fname, cachedir = Dorado.mkcacheObj(toalign, 'astrometryNet')
            path = [cachedir, fname]
            writearray = [cachedir, 'solved.fits']
            solved, wcs_header = Dorado.plate_solve(path, writearray = writearray)
            Dorado.delcacheObj( fname, 'astrometryNet')
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(wcs_header)
            Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = solved
    
    def align(self, cr, filter, alignto = None, getWCS = True, cache = False, ds = 2, ma = 5):
        '''
        align aligns a filter stack within a specified ceres object to a specified frame (default is first
        frame). align can optionally retrieve the corresponding WCS data for the aligned stack.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to align.
        filter: str
            String representation of the relevent filter.
        alignto: int
            Index of image to use as reference. Default is stack.alignto
        getWCS: boolean
            Controls whether to obtain WCS information for stack.
        cache: boolean
            Controls whether to call astrometryNet for solve or use solved result stored in cache 
            from previous run. Default is True.
        ds: int or float
            Sets the detection sigma value for astroalign.register. Default is 2.
        ma: int or float
            Sets the minimum area value for astroalign.register. Default is 5.
        '''
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        
        if alignto == None:
            alignto = series.alignTo
        else: 
            series.alignTo = alignto
        toalign = series.data[alignto]
        ## TODO :: make this use ceres.getWCS()
        if getWCS:
            if cache:
                toalign =  CCDData.read(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits', unit = Dorado.unit)
                hdulist = fits.open(Dorado.dordir / 'cache' / 'astrometryNet' / 'solved.fits') 
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(hdulist[0].header, hdulist)
                hdulist.close()
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = toalign
            else:
                fname, cachedir = Dorado.mkcacheObj(toalign, 'astrometryNet')
                path = [cachedir, fname]
                writearray = [cachedir, 'solved.fits']
                solved, wcs_header = Dorado.plate_solve(path, writearray = writearray)
                toalign = solved
                Dorado.delcacheObj( fname, 'astrometryNet')
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].wcs = WCS(wcs_header)
                Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].solved = solved
                # delete cache object
                # save solved to target

        aa_series = []
        skipped = []

        print('Aligning')
        for image in tqdm(series.data, colour = 'green'):
            # bar.update()
            try:
                img, _ = aa.register(image.data, toalign.data, detection_sigma = ds, min_area = ma)
                aaim = image
                aaim.data = img
                aa_series.append(aaim)
            except:
                skipped.append(image)
                # print('Image skipped')
        if len(skipped) != 0:
            print(len(skipped), ' images skipped.')
            ## TODO :: need to redo times and such for less ims

        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].data = aa_series
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].aligned = True
        
    def apPhot(self, cr, filter, toid, control_toid = None, shape = 21, unc = 0.1):
        '''
        apPhot performs basic aperture photometry based on photutils.aperture_photometry on a target
        within stack. Target photometry can optionally be compared to a control target within the stack 
        via the differential photometry method.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform aperture photometry on.
        filter: str
            String representation of the relevent filter.
        toid: string

        control_toid: string

        shape: int

        unc: float

        '''
        # TODO this needs a better name
        # TODO get seeing from PSF
        stack = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        # TODO if no wcs, complain alot
        w = stack.wcs
        toi = Dorado.targets[Dorado.target_keys[toid]]
        xy = w.wcs_world2pix(toi.coords.ra.deg, toi.coords.dec.deg, 1)
        ra = toi.coords.ra.deg
        dec = toi.coords.dec.deg
        pos = [(float(xy[0]), float(xy[1]))]
        aperture = CircularAperture(pos, r = shape)
        annulus_aperture = CircularAnnulus(pos, r_in = shape + 2, r_out = shape + 5)
        apers = [aperture, annulus_aperture]

        # TODO can a control be chosen? can more than one be chosen?
        if control_toid != None:
            control_toi = Dorado.targets[Dorado.target_keys[control_toid]]
            xyc = w.wcs_world2pix(control_toi.coords.ra.deg, control_toi.coords.dec.deg, 1)
            posc = [(float(xyc[0]), float(xyc[1]))]
            aperturec = CircularAperture(posc, r = shape)
            annulus_aperturec = CircularAnnulus(posc, r_in = shape + 2, r_out = shape + 5)
            apersc = [aperturec, annulus_aperturec]


        timestr = []
        exptimes = []
        ray = []
        decx = []
        x = []
        y = []
        flux = [] # NOTE this is exptime corrected
        fluxunc = []
        apsum = [] # NOTE this is raw aperture sum
        apsum_unc = []


        print('Performing photometry')
        for image in tqdm(stack.data, colour = 'green'):
            error = unc * image.data
            results = aperture_photometry(image, apers, error = error)
            bkg_mean = results['aperture_sum_1'] / annulus_aperture.area
            bkg_sum = bkg_mean * aperture.area
            results['flux_fit'] = results['aperture_sum_0'] - bkg_sum
            
            timestr.append(image.header['DATE-OBS']) # TODO This is most likely needed, but verify 
            exptimes.append(image.header['EXPTIME']) # TODO is this also needed
            ray.append(ra) # TODO is this really needed?
            decx.append(dec)
            x.append(results['xcenter'][0]) # TODO is this needed either?
            y.append(results['ycenter'][0])

            if control_toid != None:
                resultsc = aperture_photometry(image, apersc, error = error)
                bkg_meanc = resultsc['aperture_sum_1'] / annulus_aperturec.area
                bkg_sumc = bkg_meanc * aperturec.area
                resultsc['flux_fit'] = resultsc['aperture_sum_0'] - bkg_sumc

                apsum.append(results['flux_fit'][0] - resultsc['flux_fit'][0])
                flux.append((results['flux_fit'][0] - resultsc['flux_fit'][0])/image.header['EXPTIME'])
            else:
                apsum.append(results['flux_fit'][0])
                flux.append(results['flux_fit'][0]/image.header['EXPTIME'])

            fluxunc.append(1) ## TODO:: modify this to account for exposure time and control
            apsum_unc.append(1)
        times = Time(timestr)
        ts = timeSeries(times = times, flux = flux, exptimes = exptimes, x = x, y = y, 
        ra = ray, dec = decx, flux_unc = fluxunc, apsum = apsum, apsum_unc = apsum_unc)

        Dorado.targets[Dorado.target_keys[toid]].filters[filter] = len(toi.ts)
        Dorado.targets[Dorado.target_keys[toid]].ts.append(ts) 
        # TODO accomodate targets embedded in core (list of targets)
        # TODO the name for this function needs updating
        
    def mkBase(self, cr, filter, sigClip = False, minmax = False):
        '''
        mkBase stacks filter data within a ceres object to produce a base or 'average' stacked image.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to create a base stacked image for.
        filter: str
            String representation of the relevent filter.
        sigClip: boolean

        minmax: boolean

        '''
        ## TODO :: add the option to change the combination method. Right now default is 
        # sigma clipped median combination.
        series = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]]
        # toalign = series.data[series.alignTo]

        c = ccdprocx.Combiner(series.data)
        if minmax:
            c.minmax_clipping(min_clip = 0.1)
        if sigClip:
            c.sigma_clipping()
        base = c.median_combine()
        base.header['stacked'] = True
        base.header['numsubs'] = len(series.data)
        base.header['DATE-OBS'] = series.data[0].header['DATE-OBS']
        base.header['filter'] = series.filter

        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].base = base
        ## TODO :: sort out what is in the header of this base file.

        
    def calBase(self, cr, filter):
        '''
        calBase calibrates a base image for a stack by computing a 2D background for the image
        and removing it from the base image.
        
        Parameters
        ----------
        cr: string
            The relevent name string of Ceres instance to perform base image calibration on.
        filter: str
            String representation of the relevent filter.
        '''
        # TODO this needs hella optimization and direction
        # I am 99% sure this is to remove the background gradient from the stacked base image
        img = Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].base
        norm = ImageNormalize(stretch=SqrtStretch())
        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        bkg = Background2D(img, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        base = img
        base.data = img.data / bkg.background.value
        base.data[np.isnan(base.data)] = 0
        Dorado.ceres[Dorado.ceres_keys[cr]].data[Dorado.ceres[Dorado.ceres_keys[cr]].filters[filter]].base = base
        








## TODO :: figure out how to handle aligning and processing planetary images and image with low star counts







# class tessPhot:
    # def __init__(self):
    #     self.temp = None

    # def dorphot(self, cere, filter, aperture_shape, annulus_shape):
    #     stack = cere.data[cere.filters[filter]]
        
    #     if stack.ts == []:
    #         #tess_bjds, make sure to read them in as time objects for cere
    #         stack.ts = timeSeries()
    #     ts = stack.ts

    #     perture = RectangularAperture((5,4), 3,3) 
    #     annulus_aperture = RectangularAnnulus((5,4), 4, 5, 5) 
    #     # Combine the aperture
    #     aps = [aperture , annulus_aperture]

    #     for s in tqdm(range(stack.length), colour = 'green'):
    #         # Place the aperature on the current stamp and perform aperture photometry
    #         phot_table = aperture_photometry(calibrated_fluxes[s,:,:], aps)
    #         # Compute the background values from the annulus photometry
    #         bkg_sum = (phot_table['aperture_sum_1'] / annulus_aperture.area) * aperture.area
    #         # Subtract the background from main aperture
    #         phot_table['residual_aperture_sum'] = phot_table['aperture_sum_0'] - bkg_sum
    #         # Store the aperture results 
    #         aperture_sums.append(float(phot_table['aperture_sum_0']) - float(bkg_sum))
    #         # Store the flux error
    #         err.append(np.mean(flux_err[s,:,:]))
