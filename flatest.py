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

sky_flats = [CCDData.read('RED MASTER FLAT SKY.fits', unit = u.adu), CCDData.read('GREEN MASTER FLAT SKY.fits', unit = u.adu), CCDData.read('BLUE MASTER FLAT SKY.fits', unit = u.adu)]
dome_flats = [CCDData.read('RED FLAT MASTER DOME.fits', unit = u.adu), CCDData.read('GREEN MASTER FLAT DOME.fits', unit = u.adu), CCDData.read('BLUE MASTER FLAT DOME.fits', unit = u.adu)]
draco.plt_fits(np.abs(sky_flats[1].data / dome_flats[1].data), 'viridis')

# make a 2 X 3 image plot
# Rs, Gs, Bs
# Rd, Gd, Bd

## end of changes

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')