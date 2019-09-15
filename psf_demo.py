from astropy.visualization import simple_norm
from astropy.nddata import NDData
from photutils.psf import extract_stars
from photutils import EPSFBuilder
from astropy.stats import sigma_clipped_stats


import datetime
import warnings
warnings.filterwarnings('ignore')
#########

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
# from matplotlib import cm
# from astropy.visualization import astropy_mpl_style
# from astropy.utils.data import get_pkg_data_filename
# from matplotlib.colors import LogNorm
# import matplotlib.colors as colors

from astroscrappy import detect_cosmics
from astropy.nddata import CCDData
from ccdproc import Combiner
import ccdproc
import astropy.units as u
from astropy.utils.misc import isiterable
from astropy.visualization import astropy_mpl_style

import skimage.io

import os
import sys
import argparse
import astropy.io.fits as pyfits
import skimage.morphology as morph
import skimage.exposure as skie
from scipy.ndimage import median_filter

from scipy.ndimage import gaussian_filter
from skimage import data # might not be needed
from skimage import img_as_float
from skimage.morphology import reconstruction
from scipy import ndimage

from skimage.feature import blob_dog, blob_log, blob_doh
# from scipy.stats import mode

import draco
# 2018-12-04+05
# 2019-01-13+14
# 2019-08-14+15
# 2019-08-23+24
# 2019-01-04+05
directory = 'D:/draco_data/2018-12-04+05'
# star = 'BL-Cam'
star = 'YZBoo'
# star = 'DYPEG'
filename = 'FstarLum'

#header, data = draco.get_im('', filename)

# def checkdir(directory):
#     dirlist = os.listdir(directory)
#     bias = [s for s in dirlist if 'BIAS.FIT' in s]
#     flats = [s for s in dirlist if 'FLAT.FIT' in s]
#     lights = [s for s in dirlist if (np.invert('FLAT.FIT' in s)) and (np.invert('BIAS.FIT' in s)) and ('.FIT' in s) ]
#     # [s for s in dirlist if 'BIAS.FIT' in s]
#     print('\ndirlist: ', len(dirlist))
#     print('\nbias\': ', len(bias))
#     print('\nflats: ', len(flats))
#     print('\nlights: ', len(lights))
#     print('\ntotal: ', len(lights) + len(flats) + len(bias))

#     return bias, flats, lights


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

def starSeeker(data):
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
    lim = 0.85
    x2, y2 = x1[v > lim], y1[v > lim]

    return x2, y2, opt_img



def starSeeker2(data):
    mf = median_filter(data, size= 15)
    datamf = data - mf
    limg = np.arcsinh(datamf)
    limg = limg / limg.max()
    low = np.percentile(limg, 1)
    high = np.percentile(limg, 99.5)
    opt_img  = skie.exposure.rescale_intensity(limg, in_range=(low,high))


    stars =  blob_log(opt_img, max_sigma=25, min_sigma = 5, num_sigma=10, threshold=.2)
    # stars =  blob_dog(opt_img, max_sigma=30, threshold=.2)
    # Compute radii in the 3rd column.
    stars[:, 2] = stars[:, 2] * np.sqrt(2)

    y2, x2, r = stars[:, 0], stars[:, 1], stars[:, 2]

    limg = np.arcsinh(data)
    limg = limg / limg.max()
    low = np.percentile(limg, 0.1)
    high = np.percentile(limg, 99.5)
    opt_img  = skie.exposure.rescale_intensity(limg, in_range=(low,high))

    return x2, y2, r, opt_img

def plt_stars(data, x, y, r):
    plt.style.use(astropy_mpl_style)
    plt.figure()
    plt.imshow(data, cmap='viridis', vmin=0)
    up, down = plt_eye(data)
    cbar = plt.colorbar()
    cbar.set_clim(down, up)
    plt.grid(False)
    
    for i in range(len(x)):
        circlei=plt.Circle((x[i],y[i]), r[i], edgecolor='r', alpha = 0.5)
        plt.gcf().gca().add_artist(circlei)

    plt.show()

def plt_eye(data):
    mean = np.mean(data)
    std = np.std(data)
    mean_up = mean + 1.5 * std
    mean_down = mean - 1.5 * std

    return mean_up, mean_down


START_DATE_TIME = datetime.datetime.now()

print('\nStarting time: ', START_DATE_TIME)

print('\nFinding data...')
bias_list, flats_list, lights_list = draco.checkdir(directory)
#lights = draco.get_series(directory, lights_list)

print('\nStarting bias reduction')
bias, bias_headers = draco.get_series(directory, bias_list)

master_bias = draco.mast_reduc(bias)

size = master_bias.shape


print('\nStarting flat reduction')
master_i_Flat = []
flats, flat_headers = draco.get_series(directory, flats_list)
flat = draco.mast_flat(flats, master_bias)
flat = draco.norm_flat(np.array((flat.data)))

print('\nStarting light reduction')
# series, headers = draco.get_series(directory, lights_list)
series, headers = draco.get_im(directory, lights_list[0])
# series = draco.series_arith(series, '-', master_bias, series.shape[2])
# series = draco.series_arith(series, '/', flat, series.shape[2])
# sky = draco.sky_est(series)
# series = draco.series_arith(series, '-', sky, series.shape[2])
series = draco.imarith(series, '-', master_bias)
series = draco.imarith(series, '/', flat)
sky = draco.sky_est(series)
series = draco.imarith(series, '-', sky)

print('\nLight reduction complete.')


print('\nStarting star identification...')
# data = series[:, :, 10]
data = series
# size = 50
# hsize = (size - 1) / 2

# data = data[(data > hsize) & (data < (len(data[1]) -1 - hsize))]
# #data = data[(data[1] > hsize) & (data[1] < (len(data[1]) -1 - hsize))]

x, y, r, opt_img = starSeeker2(data)
# x, y, r, opt_img = starSeeker2(data)
print('\nNumber of stars: ', len(x))
print('\nStarting star field plotting...\n')
plt_stars(opt_img, x, y, r)


#########


print('\nFitting point spread functions')
# vary size based on r
size = 50
hsize = (size - 1) / 2

mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)))

stars_tbl = Table()
stars_tbl['x'] = x[mask]
stars_tbl['y'] = y[mask]

mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.)  
data = data - float(median_val)

nddata = NDData(data=data)
stars = extract_stars(nddata, stars_tbl, size=50)
nrow = 4
ncol = 4
fig, ax = plt.subplots(nrows=nrow, ncols=ncol)

ax = ax.ravel()
if (len(x) > nrow*ncol):
        ranging = nrow*ncol
else:
        ranging = len(x[mask])
for i in range(ranging):
        norm = simple_norm(stars[i], 'log', percent=99.)
        ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')

plt.show()

epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, progress_bar=True)  
epsf, fitted_stars = epsf_builder(stars)  

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))

norm = simple_norm(epsf.data, 'log', percent=99.)
plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
plt.colorbar()
plt.show()

fig, ax = plt.subplots(nrows=nrow, ncols=ncol)

ax = ax.ravel()
if (len(x) > nrow*ncol):
        ranging = nrow*ncol
else:
        ranging = len(x[mask])
for i in range(ranging):
        norm = simple_norm(fitted_stars[i], 'log', percent=99.)
        ax[i].imshow(fitted_stars[i], norm=norm, origin='lower', cmap='viridis')

plt.show()