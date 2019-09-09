from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib import cm
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from matplotlib.colors import LogNorm
import matplotlib.colors as colors

from astroscrappy import detect_cosmics
from astropy.nddata import CCDData
from ccdproc import Combiner
import ccdproc
import astropy.units as u
from astropy.utils.misc import isiterable

from scipy.stats import mode
from scipy import ndimage

import skimage.io

import os
import sys
import argparse
import skimage.morphology as morph
import skimage.exposure as skie
from scipy.ndimage import median_filter


file_path = 'C:/Users/mucep/Offline/Draco-test/'
im_prefix = 'Draco-'
c46P_prefix = '46P_'
BIAS_suffix = '.BIAS'
FLAT_prefix = 'FLAT_'
FLAT_suffix = '.FLAT'
filters = ['BLUE', 'CLEAR', 'EIGHT', 'GREEN', 'LUMINANCE', 'RED', 'SEVEN', 'SIX']

def plt_fits(image_data, cmap_str):
    plt.style.use(astropy_mpl_style)
    #print(image_data.shape)

    plt.figure()
    #plt.imshow(np.transpose(np.fliplr(image_data)), cmap='viridis')
    #plt.imshow(image_data, cmap='viridis', norm=LogNorm())
    plt.imshow(image_data, cmap=cmap_str, vmin=0)
    plt.colorbar()
    plt.grid(False)

    plt.show()
    
def plt_flat(image_data, cmap_str):
    plt.style.use(astropy_mpl_style)
    #print(image_data.shape)

    plt.figure()
    #plt.imshow(np.transpose(np.fliplr(image_data)), cmap='viridis')
    #plt.imshow(image_data, cmap='viridis', norm=LogNorm())
    plt.imshow(image_data, cmap=cmap_str, norm=LogNorm(vmin=image_data.min(), vmax=image_data.max()))
    plt.colorbar()
    plt.grid(False)

    plt.show()
    
def get_im(file_path, file_name):
    file_string = file_path + file_name + '.fits'
    # print(file_string)
    image_file = fits.open(file_string)
    # print(image_file)

    #fits.info(file_string)

    image_data = fits.getdata(file_string)

    image_header = fits.getheader(file_string, ext=0)

    return image_header, image_data

def mast_reduc(images):
    print(images.shape)
    reduc_file = []
    for i in range(images.shape[2]):
        reduc_file.append(CCDData(images[:, :, i], unit=u.adu))
    reduc_combine = Combiner(reduc_file)
    reduc_nocr = reduc_combine.average_combine()
    reduc = ccdproc.cosmicray_lacosmic(reduc_nocr)
    
    return reduc

def mast_reduc_old(file_path, im_prefix, im_suffix, im_count, cosmics):
    file_string = file_path + im_prefix+ str(1) + im_suffix + '.fits'
    temp = fits.getdata(file_string)
    temp_size = temp.shape
    #print(temp_size)

    #reduc_file = np.zeros((temp_size[0], temp_size[1], im_count))
    #reduc_file = np.zeros((im_count))
    # file_string = file_path + im_prefix+ str(1) + im_suffix + '.fits'
    # reduc_file = np.array(CCDData(fits.getdata(file_string), unit=u.adu))
    reduc_file = []
    for i in range(1, im_count):
        file_string = file_path + im_prefix+ str(i) + im_suffix + '.fits'
        #print(file_string)

        reduc_file.append(CCDData(fits.getdata(file_string), unit=u.adu))
        # if (cosmics):
        #     reduc_file[:, :, i] = detect_cosmics(fits.getdata(file_string))[1]
        # else:
        #     reduc_file[:, :, i] = fits.getdata(file_string)

    # Take average of reduc_file :

    # mean values have cosmic ray influence
    reduc_combine = Combiner(reduc_file)
    #reduc = np.median(reduc_file, axis=2)
    reduc_nocr = reduc_combine.average_combine()
    reduc = ccdproc.cosmicray_lacosmic(reduc_nocr)
    # median to remove influence of cosmic rays
    #reduc = np.median(reduc_file, axis=2)
    
    return reduc

def mast_flat_old(file_path, im_prefix, im_suffix, im_count, cosmics, master_bias):
    file_string = file_path + im_prefix+ str(1) + im_suffix + '.fits'
    temp = fits.getdata(file_string)


    reduc_file = []
    for i in range(1, im_count):
        #print('filter: ', im_prefix, ', flat: ', i)
        file_string = file_path + im_prefix+ str(i) + im_suffix + '.fits'
        data = fits.getdata(file_string)
        data = imarith(data, '-', master_bias)
        data = ndimage.median_filter(data, size=10) # , size=15
        reduc_file.append(CCDData(data, unit=u.adu))


    reduc_combine = Combiner(reduc_file)
    reduc_nocr = ccdproc.combine(reduc_file, method='average', sigma_clip=True) # , sigma_clip=True
    reduc = reduc_nocr
    print('reduced: ', im_prefix)
    plt_flat(np.array((reduc.data)), 'viridis')
    return reduc

def mast_flat(flats, master_bias):
    print(flats.shape)
    reduc_file = []
    for i in range(flats.shape[2]):
        data = flats[:, :, i]
        data = imarith(data, '-', master_bias)
        data = ndimage.median_filter(data, size = 5)
        reduc_file.append(CCDData(data, unit=u.adu))

    reduc_combine = Combiner(reduc_file)
    reduc = ccdproc.combine(reduc_file, method='average', sigma_clip=True)
    print((reduc.data).shape)
    print('reduced flat: ')
    plt_flat(reduc.data, 'viridis')
    return reduc

def stack_im(im_list):
    im = im_list[0]
    for i in range(1, len(im_list)):
        im = im + im_list[i]
    
    return im

def imarith(operand_1, operator, operand_2):
    if (isiterable(operand_1)):
        operand_1 = np.array(operand_1)
    if (isiterable(operand_2)):
        operand_2 = np.array(operand_2)
    if (operator == '+'):
        im = operand_1 + operand_2
    elif (operator == '-'):
        im = operand_1 - operand_2
    elif (operator == '/'):
        im = operand_1 / operand_2
    elif (operator == '*'):
        im = operand_1 * operand_2
    # not sure if this will work, or if needed :
    elif (operator == '^'):
        im = operand_1 ** operand_2
    
    return im

def norm_flat(flat):
    norm_factor = mode(flat, axis = None)[0]
    print("norm_factor: ", norm_factor)
    normFlat = imarith(flat, '/', norm_factor)

    return normFlat

def get_series_old(file_path, im_prefix, im_suffix, im_count, cosmics):
    file_string = file_path + im_prefix+ str(1) + im_suffix + '.fits'
    temp = fits.getdata(file_string)
    temp_size = temp.shape
    # print(temp_size)

    reduc_file = np.zeros((temp_size[0], temp_size[1], im_count))

    for i in range(1, im_count):
        file_string = file_path + im_prefix+ str(i) + im_suffix + '.fits'
        #print(file_string)
        if (cosmics):
            reduc_file[:, :, i] = detect_cosmics(fits.getdata(file_string))[1]
        else:
            reduc_file[:, :, i] = fits.getdata(file_string)

        #reduc_file[:, :, i] = detect_cosmics(fits.getdata(file_string))[1]

    return reduc_file

def series_arith(series, operator, operand, im_count):
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
        
def series_write(series, im_prefix, im_count):
    for i in range(1, im_count):
        file_nom = im_prefix + str(i) + '.fits'
        fits.writeto(file_nom, series[:,:,i])

def sky_est(im):
    sky = np.median(im)

    return sky

def checkdir(directory):
    dirlist = os.listdir(directory)
    bias = [s for s in dirlist if 'BIAS.FIT' in s]
    flats = [s for s in dirlist if 'FLAT.FIT' in s]
    lights = [s for s in dirlist if (np.invert('FLAT.FIT' in s)) and (np.invert('BIAS.FIT' in s)) and ('.FIT' in s) ]
    print('\ndirlist: ', len(dirlist))
    print('\nbias\': ', len(bias))
    print('\nflats: ', len(flats))
    print('\nlights: ', len(lights))
    print('\ntotal: ', len(lights) + len(flats) + len(bias))

    return bias, flats, lights

def get_series(directory, imlist):
    os.chdir(directory)
    file_string = imlist[0]
    temp = fits.getdata(file_string)
    temp_size = temp.shape
    # print(temp_size)

    reduc_file = np.zeros((temp_size[0], temp_size[1], len(imlist)))

    for i in range(1, len(imlist)):
        file_string = imlist[i]
        reduc_file[:, :, i] = fits.getdata(file_string)

    return reduc_file

def starSeeker(data, lim):
    # mf = median_filter(data, size=10)
    # data = draco.imarith(data, '-', mf)
    limg = np.arcsinh(data)
    limg = limg / limg.max()
    low = np.percentile(limg, 10)
    high = np.percentile(limg, 99.5)
    opt_img  = skie.exposure.rescale_intensity(limg, in_range=(low,high))
    lm = morph.local_maxima(opt_img)
    x1, y1 = np.where(lm.T == True)
    v = limg[(y1,x1)]
    x2, y2 = x1[v > lim], y1[v > lim]

    return x2, y2, opt_img

def plt_stars(data, x, y):
    plt.style.use(astropy_mpl_style)
    plt.figure()
    plt.imshow(data, cmap='viridis', vmin=0)
    plt.colorbar()
    plt.grid(False)
    
    for i in range(len(x)):
        circlei=plt.Circle((x[i],y[i]), 5, alpha = 0.5, edgecolor='r')
        plt.gcf().gca().add_artist(circlei)

    plt.show()

# function to set negative values of counts to zero

# function to remove cosmic rays on import 
# optionally import astroscrappy

# convert to RGB image utilizing (http://docs.astropy.org/en/stable/visualization/rgb.html)
# maybe even .Gif images this way (using Magneto code)

#

## End of functions ##
