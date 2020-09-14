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

__all__ = ['reduce_series', 'mastFlat', 'mastBias', 'theMask', 'sky_est', 'series_arith', 'norm_flat', 'stack_im']

def reduce_series(directory, night, imlist, flat, bias, vstar, expath = '', resize = ''):
        stardir = os.getcwd()
        imdir = directory + '/' + night + '/' + expath
        caldir = directory + night + '/wrk/calibrated/'
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
                if (resize != ''):
                                hdu.header['Resized'] = True
                                hdu.data = zoom(hdu.data, (resize, resize), order=0)
                os.chdir(caldir)
                hdu.write(vstar + expath + '-' + str(i) + '-calibrated.fit')
                series.append(hdu)
        # print(len(series))
        os.chdir(stardir)

        return series

def mastFlat(directory, night, flats, bias):
        # Allow resizing
        stardir = os.getcwd()

        path = directory + '/' + night + '/'
        os.chdir(path)

        master_flat = ccdproc.combine(flats, method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, unit=u.adu)
        master_flat.header['stacked'] = True


        caldir = directory + night + '/wrk/flats/'
        os.chdir(caldir)
        master_flat.write('master_flat.fit')

        os.chdir(stardir)

        return master_flat

def mastBias(directory, night, bias):
        # Allow resizing
        stardir = os.getcwd()

        path = directory + '/' + night + '/'
        os.chdir(path)

        master_bias = ccdproc.combine(bias, method='average', unit=u.adu)
        master_bias.meta['stacked'] = True


        caldir = directory + night + '/wrk/bias/'
        os.chdir(caldir)
        master_bias.write('master_bias.fit')

        os.chdir(stardir)

        return master_bias

def theMask(data, lx, hx, ly, hy):
        # Allow for multiple mask rectangles or circles
        mask = np.zeros(data.shape, dtype=bool)
        mask[lx:hx, ly:hy] = True
        return mask

def sky_est(im):
    # Allow for other sky definitions including masking
    sky = np.median(im)

    return sky

def series_arith(series, operator, operand, im_count):
    # mod to remove im_count and make possible to use single image
    for i in range(1, im_count):
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
    # fix imarith
    norm_factor = mode(flat, axis = None)[0]
    print("norm_factor: ", norm_factor)
    normFlat = imarith(flat, '/', norm_factor)

    return normFlat

def stack_im(im_list):
    # Is there another version of this somewhere?
    im = im_list[0]
    for i in range(1, len(im_list)):
        im = im + im_list[i]
    
    return im