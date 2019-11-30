import datetime
from astropy.time import Time
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

from scipy import interpolate
from scipy.interpolate import splprep, splev

import os
import sys


START_DATE_TIME = datetime.datetime.now()

print('\nStarting time: ', START_DATE_TIME)

print('\nFinding data...')

# directory = 'E:/draco_data/'
# night = '2019-11-19+20/aligned/fit'
# directory = 'E:/draco_data/'
# night = '2019-08-23+24/wrk/'
directory = 'E:/draco_data/2019-10-09+10/wrk/'
night = 'Transit-aligned/fitt'

os.chdir(directory)
trans_tbl = Table.read('lc_tbljune_ben.csv')
os.chdir(directory+night)

stardir = os.getcwd()
print('\n', stardir)

# trans_tbl = Table.read('lc_tbl20a.fit')

# trans_tbl = Table.read('lc_tbljune_ben.csv')

# trans_tbl = Table.read('lc_tbl15a_prv_3.csv')
# trans_tbl = Table.read('trans_tbl20.fit')
# Table.show_in_browser(trans_tbl)

xx = trans_tbl[10]
xxx = trans_tbl[6]
xxxx = trans_tbl[28]
xxxxx = trans_tbl[20]
xxxxxx = trans_tbl[3]
xxxxxxx = trans_tbl[18]
# xx = trans_tbl[26]
# xxx = trans_tbl[1]
# xxxx = trans_tbl[2]
# xxxxx = trans_tbl[0]
# xxxxxx = trans_tbl[4]
# xxxxxxx = trans_tbl[18]
# xx = trans_tbl[18]
# xxx = trans_tbl[13]
# xxxx = trans_tbl[19]
# xxxxx = trans_tbl[34]
# xxxxxx = trans_tbl[15]
# xxxxxxx = trans_tbl[12]
# xx = trans_tbl[7]
# xxx = trans_tbl[4]
# xxxx = trans_tbl[17]
# xxxxx = trans_tbl[11]
# xxxxxx = trans_tbl[12]
# xxxxxxx = trans_tbl[0]

print(xx['time-' + str(1)])
print(len(xx))
frames =115 # 351, 115, 172
x0 = np.zeros(frames)
x1 = np.zeros(frames)
x2 = np.zeros(frames)
x3 = np.zeros(frames)
x4 = np.zeros(frames)
x5 = np.zeros(frames)
ym = np.zeros(frames)
for h in range(len(x0)):
        x0[h] = xx['mag-' + str(h)]
        x1[h] = xxx['mag-' + str(h)]
        x2[h] = xxxx['mag-' + str(h)]
        x3[h] = xxxxx['mag-' + str(h)]
        x4[h] = xxxxxx['mag-' + str(h)]
        x5[h] = xxxxxxx['mag-' + str(h)]
        #datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        dt = Time(xx['time-' + str(h)], format='fits') - Time(xx['time-' + str(0)], format='fits')
        #print(dt)
        ym[h] = dt.sec
# print(x0)
# print(x1)
# print(x2)



plt.figure()
# xm = ((x0-x3) + (x0-x4)) / 2
# xm = np.log(np.abs(((x0-x1) + (x0-x2) + (x0-x3) + (x0-x4) + (x0-x5)) / 5))
xm = np.log10(np.abs(((x0-x1) + (x0-x2) + (x0-x3) + (x0-x4) + (x0-x5)) / 5))
#plt.plot(y, x0, '--', y, x1, '-', y, x2, '--')
# plt.scatter(ym, np.log(np.abs(x0-x1)), marker = 'o')
# plt.scatter(ym, np.log(np.abs(x0-x2)), marker = '+')
# # #plt.scatter(y, (x0-x3), marker = '^')
# plt.scatter(ym, np.log(np.abs(x0-x4)), marker = 'x')
# plt.scatter(ym, np.log(np.abs(x0-x5)), marker = 'D')

# spline parameters
s = 3.0 # smoothness parameter
k = 4 # spline order
nest = -1 # estimate of number of knots needed (-1 = maximal)

# find the knot points
x = xm
y = ym
# x = xm[xm > 0.515]
# y = ym[xm > 0.515]
print(len(x), ' ', len(y))
tckp,u = splprep( [x, y],s=s,k=k,nest=-1 )


# evaluate spline, including interpolated points
ynew,xnew = splev(np.linspace(0,1,400), tckp)

# f = interpolate.interp1d(y, xm)
# xnew = np.arange(0, np.max(y), 0.1)
# ynew = f(xnew)   # use interpolation function returned by `interp1d`
plt.plot(xm, ym, '.')
# plt.plot(xnew, ynew, '--')

# transit specifit limits
#plt.ylim((0.51,0.54))
plt.ylim((np.min(ym),np.max(ym)))
print(np.min(ym),np.max(ym))
plt.show()


# epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, progress_bar=True)  
# epsf, fitted_stars = epsf_builder(stars)  

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))