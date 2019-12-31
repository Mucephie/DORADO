import datetime
from astropy.time import Time
import math
import warnings
warnings.filterwarnings('ignore')
#########
import draco
import dracoOP2

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from astropy.nddata import CCDData
from ccdproc import Combiner
import ccdproc
import astropy.units as u
from astropy.utils.misc import isiterable
from astropy.visualization import astropy_mpl_style

from astropy.visualization import simple_norm
from astropy.nddata import NDData
from photutils.psf import extract_stars
from photutils import EPSFBuilder
from photutils import aperture_photometry, CircularAperture, RectangularAnnulus, RectangularAperture

import os
import sys


START_DATE_TIME = datetime.datetime.now()

print('\nStarting time: ', START_DATE_TIME)

print('\nFinding data...')

# directory = 'E:/draco_data/2017-09-20+21/wrk/aligned/'
# night = 'fits'
directory = 'E:/draco_data/2019-10-09+10/wrk/Transit-aligned/'
night = 'fitt'

bias_list, flats_list, lights_list = dracoOP2.checkdir(directory, night)

series= dracoOP2.get_series(directory, night, lights_list)


data = series[0]
print(data.header['DATE-OBS'])
t = Time(data.header['DATE-OBS'], format='fits')
x, y, r, opt_img = dracoOP2.starSeeker2(data)

print('\nNumber of stars: ', len(x))
print('\nStarting star field plotting...\n')
# dracoOP2.plt_stars(opt_img, x, y, r)


#########

# vary size based on r
size = 50
hsize = 50 

mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)))


xm = x[mask]
ym = y[mask]
rm = r[mask]
dracoOP2.plt_stars(opt_img, xm, ym, rm)


stars_tbl = Table()
stars_tbl['x'] = x[mask]
stars_tbl['y'] = y[mask]


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

trans_tbl = Table()
scatter = 20
trans_tbl['x'] = xm
trans_tbl['y'] = ym
trans_tbl['is_v'] = np.zeros(xm.shape, dtype = bool)
trans_tbl['is_c'] = np.zeros(xm.shape, dtype = bool)
mag_i = np.zeros(xm.shape)

ap_r = 20


trans_tbl['mag-0'] = mag_i


positions = (xm, ym)


for s in range(len(series)): 

        print('Frame ', str(s), ': ')
        data = series[s]

        # nddata = NDData(data=data)
        # stars = extract_stars(nddata, stars_tbl, size=50)

        for d in range(len(stars)):
                position = (xm[d], ym[d])
                aperture = RectangularAperture( position, ap_r, ap_r)
                annulus_aperture = RectangularAnnulus(position, ap_r, ap_r + 2, ap_r + 2)
                aps = [aperture , annulus_aperture]
                phot_table = aperture_photometry(data, aps)
                # Background subtract
                bkg_sum = (phot_table['aperture_sum_1'] / annulus_aperture.area) * aperture.area
                #Subtracting background from main aperture
                #print(phot_table['aperture_sum_0'] - bkg_sum)
                mag_i[d] = float((phot_table['aperture_sum_0'] - bkg_sum).value)


        trans_tbl['sum-' + str(s)] = np.zeros(trans_tbl['x'].shape)
        trans_tbl['mag-' + str(s)] = np.zeros(trans_tbl['x'].shape)

        t = []
        for p in range(len(mag_i)):
                trans_tbl['sum-' + str(s)][p] = mag_i[p]
                t.append(str(data.header['DATE-OBS']))
                trans_tbl['mag-' + str(s)][p] = -2.5 * np.log10(mag_i[p])
        trans_tbl['time-' + str(s)] = t


stardir = os.getcwd()
print('\n', stardir)

trans_tbl.write('lc_tblkelt16b.csv')

 

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))

