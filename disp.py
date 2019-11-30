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
import astroalign as aa

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
import waldo2funcs as w2f

# FOV 60cm      9.67 x 6.44 arcmin

## testing
START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)

home_dir = os.getcwd()
print('\nHome dir: ', home_dir)
target = input('Enter the target name (e.g DYPEG): ')
usetdy = input('Reduce todays data (y/n)? ')
if (usetdy == 'y'):
        night = dracoOP2.get_night()
elif (usetdy == 'z'):
        night = '2019-08-14+15'
elif (usetdy == 'x'):
        night = '2019-08-23+24'
else:
        night = input('Enter night code (YYYY-MM-DD+DD): ')

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
directory = 'D:/draco_data/'
datadir = directory + night + '/wrk/calibrated/'

bias_list, flats_list, lights_list = dracoOP2.checkdir(directory, night, '/wrk/calibrated/')

series = dracoOP2.get_series(directory, night, lights_list, '/wrk/calibrated/', unit = None)




draco.plt_fits(series[106], 'viridis')

## end of changes

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')
