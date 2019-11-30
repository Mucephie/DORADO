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

from astroscrappy import detect_cosmics
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
from photutils import aperture_photometry, CircularAperture

import os
import sys


START_DATE_TIME = datetime.datetime.now()

print('\nStarting time: ', START_DATE_TIME)

print('\nFinding data...')

# directory = 'E:/draco_data/'
# night = '2019-11-19+20/aligned/fit'
# night = '2019-08-23+24/wrk/aligned/fit'
directory = 'E:/draco_data/2019-10-09+10/wrk/Transit-aligned/'
night = 'fitt'

bias_list, flats_list, lights_list = dracoOP2.checkdir(directory, night)
#lights = draco.get_series(directory, lights_list)
series= dracoOP2.get_series(directory, night, lights_list)


data = series[0]
print(data.header['DATE-OBS'])
t = Time(data.header['DATE-OBS'], format='fits')
x, y, r, opt_img = dracoOP2.starSeeker2(data)
# x, y, r, opt_img = starSeeker2(data)
print('\nNumber of stars: ', len(x))
print('\nStarting star field plotting...\n')
# dracoOP2.plt_stars(opt_img, x, y, r)


#########


print('\nFitting point spread functions')
# vary size based on r
size = 50
hsize = 50 # (size) / 2

# mask = dracoOP2.theMask(data, int(hsize), int((data.shape)[0]-hsize), int(hsize), int((data.shape)[1]-hsize)) # ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)))
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
trans_tbl['is_v'] = np.zeros(xm.shape, dtype=bool)
trans_tbl['is_c'] = np.zeros(xm.shape, dtype=bool)
mag_i = np.zeros(xm.shape)

aperture_radius =15.0

for s in range(len(stars)):
        position = (xm[s], ym[s])
        aperture = CircularAperture(position, r = aperture_radius)
        phot_table = aperture_photometry(data, aperture)
        # mag_i[s] = np.sum(stars[s])

trans_tbl['mag-0'] = mag_i

#print(data.)
#aperture_radius = 7.0
# positions = (results['x'], results['y'])
# apertures = CircularAperture(positions, r=aperture_radius)
# phot_table = aperture_photometry(data, apertures)
# results['aperture_sum'] = phot_table['aperture_sum']
# # add a col with calculation for instrumental mag

# results['instrumental_mag'] = results.apply(lambda x: -2.5 * math.log10(x['aperture_sum']), axis = 1)

positions = (xm, ym)


for s in range(len(series)): # len(series)

        print('Frame ', str(s), ': ')
        data = series[s]
        # x, y, r, opt_img = dracoOP2.starSeeker2(data)
        #print('\nNumber of stars: ', len(x))
        #mask = ((x > hsize) & (x < (data.shape[1] - hsize)) & (y > hsize) & (y < (data.shape[0] - hsize)))
        #xm = x[mask]
        #print('xm size: ', xm.size)
        #print('data size: ', data.shape)
        #ym = y[mask]
        #rm = r[mask]
        nddata = NDData(data=data)
        stars = extract_stars(nddata, stars_tbl, size=50)

        # positions = (x, y)
        apertures = CircularAperture(positions, r = aperture_radius)
        phot_table = aperture_photometry(data, apertures)
        mag_i = phot_table['aperture_sum']
        # i_mag = -2.5 * math.log10(mag_i)

        # mag = np.zeros(xm.shape)
        # for q in range(len(stars)):
        #         mag[q] = np.sum(stars[q])
        trans_tbl['sum-' + str(s)] = np.zeros(trans_tbl['x'].shape)
        trans_tbl['mag-' + str(s)] = np.zeros(trans_tbl['x'].shape)
        # trans_tbl['time-' + str(s)] = np.array(trans_tbl['x'].shape, dtype='')
        t = []
        for p in range(len(mag_i)):
                trans_tbl['sum-' + str(s)][p] = mag_i[p].value
                t.append(str(data.header['DATE-OBS']))
                if (mag_i[p] > 0):
                        trans_tbl['mag-' + str(s)][p] = -2.5 * np.log10(mag_i[p].value)
        trans_tbl['time-' + str(s)] = t

        
        # for w in range(len(trans_tbl['x'])):
        #         count = 0
        #         for c in range(len(stars)):
        #                 if (np.abs(x[c]-trans_tbl['x'][w]) < scatter) and (np.abs(y[c]-trans_tbl['y'][w]) < scatter) :
        #                         count = count + 1
        #                         trans_tbl['sum-' + str(s)][w] = mag_i[c].value
        #                         trans_tbl['mag-' + str(s)][w] = -2.5 * math.log10(mag_i[c].value)
        #         if count > 1:
        #                 print('WARNING :: More than one star found for star-', str(w), ' in frame ', str(s), ' !!')
        #                 print(count)
        #         if count < 1:
        #                 print('WARNING :: Less than one star found for star-', str(w), ' in frame ', str(s), ' !!')

stardir = os.getcwd()
print('\n', stardir)

trans_tbl.write('lc_tbljune_ben.csv')

# plt.figure()
# plt.plot(trans_tbl[2, 4:], np.linspace(0,len(trans_tbl[2, 4:])))
# plt.show()


# epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, progress_bar=True)  
# epsf, fitted_stars = epsf_builder(stars)  

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))

# norm = simple_norm(epsf.data, 'log', percent=99.)
# plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
# plt.colorbar()
# plt.show()

# fig, ax = plt.subplots(nrows=nrow, ncols=ncol)

# ax = ax.ravel()
# if (len(x) > nrow*ncol):
#         ranging = nrow*ncol
# else:
#         ranging = len(x[mask])
# for i in range(ranging):
#         norm = simple_norm(fitted_stars[i], 'log', percent=99.)
#         ax[i].imshow(fitted_stars[i], norm=norm, origin='lower', cmap='viridis')

# plt.show()