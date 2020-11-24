import warnings
warnings.filterwarnings('ignore')

## file imports
import os
from astropy.io import fits 
from astropy.nddata import CCDData
from astropy.table import Table
from astropy.stats import mad_std

## Processing imports
import numpy as np
import ccdproc
from scipy.ndimage import zoom

## photometry imports
from astropy import units as u

# __all__ = ['reduce_series', 'mastFlat', 'mastBias', 'theMask', 'sky_est', 'series_arith', 'norm_flat', 'stack_im']

def reduce_series(imdir, imlist, flat, bias, target, resize = '', caldir = ''):
        # directory, night, imlist, flat, bias, vstar, expath = '', resize = ''
        """
        plate_solve takes fits image file data and the corresponding file string to the data and calls nova.astrometry.net to obtain and then integrate WCS into the HDU.

        Parameters
        ----------
        imdir: filestring
                Path to the image data directory.
        caldir: filestring
                Path to store the calibrated image data.
        imlist: array[string]
                array of filestrings of data to reduce.
        flat: CCDdata
                CCDdata of the flatfield.
        bias: CCDdata
                CCDdata of the bias.
        target: string
                Name of target to use in file output name.
        resize: float
                scale factor to resize image data. default is None
        caldir: filestring
                Path to store the calibrated image data.
            
        Returns
        -------
        series: array[CCDdata]
                The CCDdata array containing the reduced images.
        """
        # add ability to use series already loaded into memory.
        stardir = os.getcwd()
        # imdir = directory + '/' + night + '/' + expath
        # caldir = directory + night + '/wrk/calibrated/'
        if caldir == '':
                caldir = imdir + '/wrk/calibrated/'
        os.chdir(caldir)
        series = []
        for i in range(len(imlist)):
                os.chdir(imdir)
                hdu = CCDData.read(imlist[i], unit=u.adu)
                hdu = ccdproc.ccd_process(hdu, master_bias = bias, master_flat = flat)
                hdu = ccdproc.cosmicray_lacosmic(hdu, sigclip=5)
                hdu.header['bias corrected'] = True
                hdu.header['flat corrected'] = True
                hdu.header['cosmicray corrected'] = True
                hdu.header['calibrated'] = True
                # may want to move this to before reduction step
                if (resize != ''):
                                hdu.header['Resized'] = True
                                hdu.data = zoom(hdu.data, (resize, resize), order=0)
                os.chdir(caldir)
                # mod to use imlist imlist[i]
                hdu.write(target + '-' + str(i) + '_c.fit', overwrite = True)
                series.append(hdu)
        # print(len(series))
        os.chdir(stardir)

        return series

def mastFlat(directory, flats, bias, caldir = ''):
        """
        mastFlat takes a datadirectory string, a list of flats and a master bias image to construct 
        a master flatfield image.

        Parameters
        ----------
        directory: filestring
                Path to the image data directory.
        flats: array[CCDdata]
                array of raw flatfields.
        bias: CCDdata
                CCDdata of the bias.
        caldir: filestring
                Path to store the calibrated image data.

        Returns
        -------
        master_flat: CCDdata
                The combined master flatfield image.
        """
        # give option to median stack
        # give option to specify a filter for output.
        # directory, night, flats, bias
        # Allow resizing
        stardir = os.getcwd()

        path = directory # + '/' + night + '/'
        os.chdir(path)

        master_flat = ccdproc.combine(flats, method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, unit=u.adu)
        master_flat.header['stacked'] = True

        if caldir == '':
                caldir = directory + '/wrk/flats/'
        os.chdir(caldir)
        master_flat.write('master_flat.fit', overwrite = True)

        os.chdir(stardir)

        return master_flat

def mastBias(directory, bias, caldir = ''):
        """
        mastBias takes a datadirectory string, a list of bias images to construct 
        a master bias image.

        Parameters
        ----------
        directory: filestring
                Path to the image data directory.
        bias: array[CCDdata]
                array of raw bias images.
        caldir: filestring
                Path to store the calibrated image data.

        Returns
        -------
        master_bias: CCDdata
                The combined master bias image.
        """
        # directory, night, bias
        # Allow resizing
        # Allow specification of median or mean
        stardir = os.getcwd()

        path = directory # + '/' + night + '/'
        os.chdir(path)

        master_bias = ccdproc.combine(bias, method='average', unit=u.adu)
        master_bias.meta['stacked'] = True

        if caldir == '':
                caldir = directory +  '/wrk/bias/'
        os.chdir(caldir)
        master_bias.write('master_bias.fit', overwrite = True)

        os.chdir(stardir)

        return master_bias

def theMask(data, lx, hx, ly, hy):
        """
        theMask takes an image and bound parameters to mask the outter edges of 
        the image in the case of imperfections. Data suffering from descecant 
        failure can be made usable via this.

        Parameters
        ----------
        data: CCDdata
                CCDdata of the image to be masked.
        lx, hx, ly, hy: integers
                The low and high x coordinates and low and high y coordinates
                respectively.
        Returns
        -------
        mask: array[Boolean]
                The masked array for application to the data image in further functions.
                
        """
        # Allow for multiple mask rectangles or circles
        # Star masking
        mask = np.zeros(data.shape, dtype=bool)
        mask[lx:hx, ly:hy] = True
        return mask

def sky_est(im):
        """
        sky_est makes a crude estimate of an images sky value based on a median
        method.

        Parameters
        ----------
        im: CCDdata
                CCDdata of the image to be solved.

        Returns
        -------
        sky: float
                image sky value.
        """
        # Allow for other sky definitions including masking
        # Use one from Proj1
        sky = np.median(im)

        return sky

def series_arith(series, operator, operand):
        """
        series_arith allows one to perform mathematical(algebraic) on a series of images.

        Parameters
        ----------
        series: array[CCDdata]
                Array of CCDdata of the image(s) to perform operations on.
        operator: string
                string of the algebaric operator.
        operand: CCDdata
                image to be used as an operand.

        Returns
        -------
        series: array[CCDdata]
                Array of CCDdata .
        """
        # mod to check datatype using type()
        # mod to remove im_count and make possible to use single image
        # mod to accomodate CCDdata object
        for i in range(len(series[0,0,:])):
            if (operator == '+'):
                series[:, :, i] = series[:, :, i]  + operand
            elif (operator == '-'):
                series[:, :, i] = series[:, :, i] - operand
            elif (operator == '/'):
                series[:, :, i] = series[:, :, i]  / operand
            elif (operator == '*'):
                series[:, :, i] = series[:, :, i]  * operand
        
        return series

def norm_flat(flat):
        """
            norm_flat takes a flatfield image and normalizes it for use in reduction.

            Parameters
            ----------
            flat: CCDdata
                CCDdata of the flatfield.

            Returns
            -------
            normFlat: CCDdata
                    The CCDdata object of the normalized flat.
        """
        # fix series_arith
        # beta test this function
        norm_factor = np.mode(flat, axis = None)[0]
        print("norm_factor: ", norm_factor)
        normFlat = series_arith([flat], '/', norm_factor)

        return normFlat

def stack_im(im_list):
        """
            stack_im stacks a series of images into a single frame.

            Parameters
            ----------
            imlist: array[string]
                array of filestrings of data to stack.

            Returns
            -------
            im: CCDdata
                    The stacked CCDdata image.
        """
        # Is there another version of this somewhere?
        # update for CCDdata
        # Add different methods for stacking (im combine)
        im = im_list[0]
        for i in range(1, len(im_list)):
            im = im + im_list[i]
        
        return im