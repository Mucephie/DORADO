import warnings
warnings.filterwarnings('ignore')

## file imports
import os
from astropy.io import fits 
from astropy.nddata import CCDData

## Processing imports
import numpy as np
import ccdproc
from astropy.stats import mad_std

## photometry imports


## plotting imports
import matplotlib.pyplot as plt
import matplotlib as mpl

## other imports
import datetime




## functions
def mkwrkdir(datadir, night):
        stardir = os.getcwd()

        path = datadir + '/' + night
        os.chdir(path)
        os.mkdir(path + '/wrk')
        path = path + '/wrk'
        os.mkdir(path + '/bias')
        os.mkdir(path + '/flats')
        os.mkdir(path + '/lights')
        os.mkdir(path + '/calibrated')
        os.mkdir(path + '/stars')
        os.mkdir(path + '/results')

        os.chdir(stardir)

def get_night():
      # currently does not support first/last of the month
      date = datetime.date.today()
      year = date.year
      month = date.month
      date2 = date.day
      date1 = date2 - 1
      night = str(year) + '-' + str(month) + '-' + str(date1) + '+' + str(date2)
      return night

def checkdir(directory, night):
        path = directory + '/' + night + '/'
        dirlist = os.listdir(path)
        bias = [s for s in dirlist if 'BIAS.FIT' in s]
        if len(bias)==0:
                bias = [s for s in dirlist if 'Bias.fit' in s]

        flats = [s for s in dirlist if 'FLAT.FIT' in s]
        if len(flats)==0:
                flats = [s for s in dirlist if 'FlatField.fit' in s]
        lights = [s for s in dirlist if (np.invert('FLAT.FIT' in s)) and (np.invert('BIAS.FIT' in s)) and ('.FIT' in s) ]
        if len(lights)==0:
                lights = [s for s in dirlist if (np.invert('FlatField.fit' in s)) and (np.invert('Bias.fit' in s)) and ('.fit' in s) ]
        print('\ndirlist: ', len(dirlist))
        print('\nbias\': ', len(bias))
        print('\nflats: ', len(flats))
        print('\nlights: ', len(lights))
        print('\ntotal: ', len(lights) + len(flats) + len(bias))

        return bias, flats, lights

def get_series(directory, night, imlist):
        stardir = os.getcwd()

        path = directory + '/' + night + '/'
        os.chdir(path)
        series = ccdproc.ImageFileCollection(filenames=imlist)

        os.chdir(stardir)

        return series

def write_series(directory, night, series, vstar):
        stardir = os.getcwd()
        caldir = directory + night + '/calibrated'
        os.chdir(caldir)
        for hdu, i in series.hdus():
                hdu.header['calibrated'] = True
                hdu.writeto(vstar + str(i) + '.cal')
        os.chdir(stardir)

def mastBias(directory, night, bias):
        master_bias = ccdproc.combine(bias, method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
        master_bias.meta['stacked'] = True

        stardir = os.getcwd()

        caldir = directory + night + '/bias'
        os.chdir(caldir)
        master_bias.write('master_bias.fit')

        os.chdir(stardir)

        return master_bias

def mastFlat(directory, night, flats, bias):
        for hdu, i in flats.hdus():
                hdu = ccdproc.subtract_bias(hdu, bias)
                hdu.header['bias corrected'] = True
                hdu.header['calibrated'] = True
                

        master_flat = ccdproc.combine(flats, method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
        master_flat.meta['stacked'] = True

        stardir = os.getcwd()

        caldir = directory + night + '/flats'
        os.chdir(caldir)
        master_flat.write('master_flat.fit')

        os.chdir(stardir)

        return master_flat

def reduce_series(directory, night, series, flats, bias):
        for hdu, i in series.hdus():
                        hdu = ccdproc.subtract_bias(hdu, bias)
                        hdu.header['bias corrected'] = True
                        hdu.header['flat corrected'] = True
                        hdu.header['calibrated'] = True

## testing
START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)

home_dir = os.getcwd()
print('\nHome dir: ', home_dir)

night = get_night()
directory = 'C:/Users/mucep/Offline/Draco-test/'
mkwrkdir(directory, night)

bias_list, flats_list, lights_list = checkdir(directory)

biaseries = get_series(directory, night, bias_list)
bias = mastBias(directory, night, biaseries)

flats = get_series(directory, night, flats_list, bias)
flat = mastFlat(directory, night, flats, bias)

#ccdproc.ccd_process(

