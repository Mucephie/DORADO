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

target = 'null'

night = '2019-10-28+29/June/Kelt-16b/wrk/Master'

directory = 'D:/draco_data'


path = directory + '/' + night + '/'
os.chdir(path)

bias = CCDData.read('mBIAS.fit', unit=u.adu)
print(bias.data.shape)
flat_R = CCDData.read('mRFlat.fit', unit=u.adu)
flat_V = CCDData.read('mVFlat.fit', unit=u.adu)
flat_B = CCDData.read('mBFlat.fit', unit=u.adu)

bias.data = zoom(bias.data, (2, 2), order=0)

scanrow = np.zeros((1, bias.data.shape[1]))
bias.data = np.vstack((bias.data,scanrow)) # = biastemp.data

flat_R.data = zoom(flat_R.data, (2, 2), order=0)
flat_R.data = np.vstack((flat_R.data,scanrow))
flat_B.data = zoom(flat_B.data, (2, 2), order=0)
flat_B.data = np.vstack((flat_B.data,scanrow))
flat_V.data = zoom(flat_V.data, (2, 2), order=0)
flat_V.data = np.vstack((flat_V.data,scanrow))


fits.writeto('mBias1.fit', bias.data)
fits.writeto('mBFlat1.fit', flat_B.data)
fits.writeto('mVFlat1.fit', flat_V.data)
fits.writeto('mRFlat1.fit', flat_R.data)