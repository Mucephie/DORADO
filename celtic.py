import warnings
warnings.filterwarnings('ignore')

## file imports
import os
from astropy.io import fits 
from astropy.nddata import CCDData
from astropy.table import Table
from astropy.nddata import NDData
# from pathlib import Path

## Processing imports
import numpy as np
import ccdproc
from scipy.ndimage import zoom

## photometry imports
from astropy import units as u


## plotting imports
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.visualization import astropy_mpl_style
from matplotlib import cm

## other imports
import datetime
import dracoOP2
import draco
import time

START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)

home_dir = os.getcwd()
print('\nHome dir: ', home_dir)

# target = input('Enter the target name (e.g DYPEG): ')
target = 'NGC7640'
# usetdy = input('Reduce todays data (y/n)? ')
# if (usetdy == 'y'):
#         night = dracoOP2.get_night()
# elif (usetdy == 'z'):
#         night = '2019-10-07+08'
# elif (usetdy == 'x'):
#         night = '2019-10-09+10'
# elif (usetdy == 'j'):
#         night = 'June'
# else:
#         night = input('Enter night code (YYYY-MM-DD+DD): ')
night = 'June'
# # data directory: D:/draco_data/
# # test directory: C:/Users/mucep/Offline/Draco-test/
# ## nights
# # 2018-10-9+10 GALAXIES
# # 2018-12-04+05 DYPEG $$
# # 2018-12-18+19 46P & FSTAR
# # 2019-01-04+05 BLCAM $$
# # 2019-01-13+14 BLCAM $$
# # 2019-08-14+15 YZBOO $$
# # 2019-08-23+24 XXCYG $$
directory = 'D:/Celtic/'


dracoOP2.mkwrkdir(directory, night)

bias_list, flats_list, lights_list = dracoOP2.checkdir(directory, night)

# biaseries = dracoOP2.get_series(directory, night, bias_list)
# bias = dracoOP2.mastBias(directory, night, bias_list)


flats = dracoOP2.get_series(directory, night, flats_list)
path = directory + '/' + night + '/'
os.chdir(path)

bias = CCDData.read('BIAS.fit', unit=u.adu)
flat_R = CCDData.read('R_MFLAT.fit', unit=u.adu)
flat_V = CCDData.read('V_MFLAT.fit', unit=u.adu)
flat_B = CCDData.read('B_MFLAT.fit', unit=u.adu)
# flat = dracoOP2.mastFlat(directory, night, flats_list, bias)

# might need to do a .data addition
bias3 = CCDData.read('BIAS.fit', unit=u.adu)
bias3.data = zoom(bias3.data, (1/1.5, 1/1.5), order=0)
flat_R3 = CCDData.read('R_MFLAT.fit', unit=u.adu)
flat_R3.data = zoom(flat_R3.data, (1/1.5, 1/1.5), order=0)
flat_B3 = CCDData.read('V_MFLAT.fit', unit=u.adu)
flat_B3.data = zoom(flat_B3.data, (1/1.5, 1/1.5), order=0)
flat_V3 = CCDData.read('B_MFLAT.fit', unit=u.adu)
flat_V3.data = zoom(flat_V3.data, (1/1.5, 1/1.5), order=0)


print(bias.data.shape)
print(bias3.data.shape)

dumb_list, dumber_list, R_list = dracoOP2.checkdir(directory, night, 'R')
dracoOP2.mkwrkdir(directory, night, 'R')
R_series = dracoOP2.get_series(directory, night, R_list,'R')
print(R_series[0].data.shape)
R_series = dracoOP2.reduce_seriesR(directory, night, R_list, flat_R3, bias3, target, 1.5, 'R')


dumb_list, dumber_list, B_list = dracoOP2.checkdir(directory, night, 'B')
dracoOP2.mkwrkdir(directory, night, 'B')
B_series = dracoOP2.get_series(directory, night, B_list, 'B')
print(B_series[0].data.shape)
B_series = dracoOP2.reduce_series(directory, night, B_list, flat_B, bias, target, 'B')

dumb_list, dumber_list, V_list = dracoOP2.checkdir(directory, night, 'V')
dracoOP2.mkwrkdir(directory, night, 'V')
V_series = dracoOP2.get_series(directory, night, V_list, 'V')
print(V_series[0].data.shape)
V_series = dracoOP2.reduce_series(directory, night, V_list, flat_V, bias, target, 'V')


# make a 2 X 3 image plot
# Rs, Gs, Bs
# Rd, Gd, Bd

## end of changes

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')