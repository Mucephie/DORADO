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
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from astropy.table import Table
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

import astropy.units as u

import lightkurve as lk
from lightkurve.lightcurve import TessLightCurve as tlc

import os
import sys


START_DATE_TIME = datetime.datetime.now()

print('\nStarting time: ', START_DATE_TIME)

print('\nFinding data...')

# directory = 'E:/draco_data/'
# night = '2019-11-19+20/aligned/fit'
# directory = 'E:/draco_data/'
# night = '2019-08-23+24/wrk/'
directory = 'D:/draco_data/2017-09-20+21/wrk/'
night = 'aligned/fits'

# directory = 'E:/draco_data/2019-10-09+10/wrk/'
# night = 'Transit-aligned/fitt'

os.chdir(directory)
# trans_tbl = Table.read('lc_tblkelt16b.csv')
tbl = Table.read('lc_tblxxcyg.csv')
os.chdir(directory+night)

stardir = os.getcwd()
print('\n', stardir)

# xxcyg
tid = 2236404589312109440 # GAIA DR2
xx = tbl[11]
xxx = tbl[1]
xxxx = tbl[2]
xxxxx = tbl[17]
xxxxxx = tbl[6]
xxxxxxx = tbl[20]
# kelt-16b
# xx = trans_tbl[6]
# xxx = trans_tbl[10]
# xxxx = trans_tbl[1]
# xxxxx = trans_tbl[3]
# xxxxxx = trans_tbl[7]
# xxxxxxx = trans_tbl[11]

print(xx['time-' + str(1)])
print(len(xx))
frames = 160 # 115 160
x0 = np.zeros(frames)
x1 = np.zeros(frames)
x2 = np.zeros(frames)
x3 = np.zeros(frames)
x4 = np.zeros(frames)
x5 = np.zeros(frames)
ym = np.zeros(frames)
dt = []
for h in range(len(x0)):
        x0[h] = xx['sum-' + str(h)]
        x1[h] = xxx['sum-' + str(h)]
        x2[h] = xxxx['sum-' + str(h)]
        x3[h] = xxxxx['sum-' + str(h)]
        x4[h] = xxxxxx['sum-' + str(h)]
        x5[h] = xxxxxxx['sum-' + str(h)]
        #datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        dt.append(xx['time-' + str(h)]) # - Time(xx['time-' + str(0)], format='fits') '%Y-%m-%dT%H:%M:%S.%f'
adt = dt
dt = Time(dt, format = 'fits')
# ts = Time(dt)

# for ti in range(len(dt)):
#         dt[ti] = dt[ti] - dt[0]
fig, ax = plt.subplots()
# xm = ((x0-x3) + (x0-x4)) / 2
# xm = np.log(np.abs(((x0-x1) + (x0-x2) + (x0-x3) + (x0-x4) + (x0-x5)) / 5))
flux = ((x0/x1) + (x0/x2) + (x0/x3) + (x0/x4) + (x0/x5)) / 5
#plt.plot(y, x0, '--', y, x1, '-', y, x2, '--')
plt.scatter(dt.value, (x0/x1)/np.linalg.norm((x0/x1)), marker = 'o')
plt.scatter(dt.value,  (x0/x2)/np.linalg.norm((x0/x2)), marker = '+')
plt.scatter(dt.value,  (x0/x3)/np.linalg.norm((x0/x3)), marker = '^')
plt.scatter(dt.value,  (x0/x4)/np.linalg.norm((x0/x4)), marker = 'x')
plt.scatter(dt.value,  (x0/x5)/np.linalg.norm((x0/x5)), marker = 'D')
# print(np.max(x0))
# print(np.max(x1))
# print(np.max(x2))
# print(np.max(x3))
# print(np.max(x4))
# print(np.max(x5))





N = len(flux)
# print(dt)
xt = np.arange(0, len(flux), 1)
# ax.plot(dt.value, flux, '.')
plt.xlabel('Time')
plt.ylabel('Mean Differential Flux Ratio')
fig.suptitle('xxcyg - 2017-09-20+21')
ax.xaxis.set_major_locator(MultipleLocator(20))
# For the minor ticks, use no labels; default NullFormatter.
ax.xaxis.set_minor_locator(MultipleLocator(5))
plt.grid

plt.gcf().autofmt_xdate()  # orient date labels at a slant 

# plt.xlim(Time('2019-10-09T01:51:10.890').to_datetime(), Time('2019-10-09T01:51:10.890').to_datetime())
# transit specifit limits
# plt.ylim((0, 0.17))
#plt.ylim((np.min(ym),np.max(ym)))
# print(np.min(flux), np.max(flux))
plt.show()

fluid_flux =  (x0/x2)/np.linalg.norm((x0/x2))
# lc = lk.LightCurve(time = dt.jd, flux = flux, flux_err = None, flux_unit = None, time_format = None, time_scale = 'tdb',  targetid = tid, label = 'XXCYGNI')
lc = lk.LightCurve(time = dt.jd, flux = fluid_flux, flux_err = None, flux_unit = None, time_format = None, time_scale = 'tdb',  targetid = tid, label = 'XXCYGNI')
# Print the new lightcurve object properties for record
# lc.show_properties()
# Clean the series of outliers > 7 sigma from the series
lc = lc.remove_outliers(sigma=2.8)
# Remove the downward flux trend expirienced by Tess throught an orbit
# lc = lc.flatten(window_length=1001)
lc.scatter()
plt.show()

lc_periodogram = lc.to_periodogram()
lc_periodogram.plot()

print('Max power : ', lc_periodogram.max_power)
print('Frequency at max power : ', lc_periodogram.frequency_at_max_power)
print('Period at max power : ', lc_periodogram.period_at_max_power)
# print(lc_periodogram.depth_at_max_power)
# print(lc_periodogram.duration_at_max_power)
# print(lc_periodogram.transit_time_at_max_power)
# print(lc_periodogram.period)

plt.show()
END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))