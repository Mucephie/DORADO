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

    fits.info(file_string)

    image_data = fits.getdata(file_string)

    image_header = fits.getheader(file_string, ext=0)

    return image_header, image_data

def mast_reduc(file_path, im_prefix, im_suffix, im_count, cosmics):
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

def mast_flat(file_path, im_prefix, im_suffix, im_count, cosmics, master_bias):
    file_string = file_path + im_prefix+ str(1) + im_suffix + '.fits'
    temp = fits.getdata(file_string)
    # temp_size = temp.shape

    reduc_file = []
    for i in range(1, im_count):
        print('filter: ', im_prefix, ', flat: ', i)
        file_string = file_path + im_prefix+ str(i) + im_suffix + '.fits'
        data = fits.getdata(file_string)
        data = imarith(data, '-', master_bias)
        data = ndimage.median_filter(data, size=5) # , size=15
        # reduc_file.append(data)
        #plt_flat(data, 'viridis')
        reduc_file.append(CCDData(data, unit=u.adu))

    # reduc = reduc / mean
    reduc_combine = Combiner(reduc_file)
    reduc_nocr = ccdproc.combine(reduc_file, method='median') # , sigma_clip=True
    # reduc = ccdproc.cosmicray_lacosmic(reduc_nocr)
    reduc = reduc_nocr
    print('reduced: ', im_prefix)
    plt_flat(np.array((reduc.data)), 'viridis')
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
    norm_factor = mode(flat)[0]
    norm_factor = mode(norm_factor, 0)[0]
    normFlat = imarith(flat, '/', norm_factor)
    # normFlat = CCDData(normFlat, unit = u.adu)
    # normFlat = ccdproc.cosmicray_lacosmic(normFlat)
    # normFlat = np.asarray(normFlat)

    return normFlat

def get_series(file_path, im_prefix, im_suffix, im_count, cosmics):
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

# function to set negative values of counts to zero

# function to remove cosmic rays on import 
# optionally import astroscrappy

# convert to RGB image utilizing (http://docs.astropy.org/en/stable/visualization/rgb.html)
# maybe even .Gif images this way (using Magneto code)

#

## End of functions ##
