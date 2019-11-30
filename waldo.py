import warnings
warnings.filterwarnings('ignore')

## file imports
import os
from astropy.io import fits 
from astropy.nddata import CCDData
from astropy.table import Table
# from pathlib import Path

## Processing imports
import numpy as np
import ccdproc
from astropy.stats import mad_std
import astroalign

## photometry imports
from astropy import units as u

## plotting imports
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm

## other imports
import datetime
import dracoOP2



## testing
START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)

home_dir = os.getcwd()
print('\nHome dir: ', home_dir)
target = input('Enter the target name (e.g DYPEG): ')
usetdy = input('Reduce todays data (y/n)? ')
if (usetdy == 'y'):
        night = dracoOP2.get_night()
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
directory = 'C:/Users/mucep/Offline/Draco-test/'
datadir = directory + night + '/wrk/calibrated/'

bias_list, flats_list, lights_list = dracoOP2.checkdir(directory, night, '/wrk/calibrated/')

series = dracoOP2.get_series(directory, night, lights_list, '/wrk/calibrated/', unit = None)

# x, y, r, opt_img = dracoOP2.starSeeker2(series[20])

# dracoOP2.plt_stars(series[20], x, y, r)

colors = ['r', 'g', 'b', 'y', 'cyan', 'w', 'm']

stardir = directory + night + '/wrk/stars/'
os.chdir(stardir)

# field = Table.read(str(target) + '-separations.fit')
# fieldinfo = Table.read(str(target) + '-fieldinfo.fit')

# IDd, isStar = dracoOP2.starange(x, y, r, field, fieldinfo, 20)
# # IDd.show_in_browser()
# xis = []
# yis = []
# ris = []
# for d in range(len(isStar)):
#         xis.append(x[isStar[d]])
#         yis.append(y[isStar[d]])
#         ris.append(r[isStar[d]])
# dracoOP2.plt_stars(series[20], xis, yis, ris)

## try tw0

p, (pos_0, pos_1) = astroalign.find_transform(series[0].data, series[20].data)

for (x1, y1), (x2, y2) in zip(pos_0, pos_1):
    print("({:.2f}, {:.2f}) in source --> ({:.2f}, {:.2f}) in target"
          .format(x1, y1, x2, y2))


END_DATE_TIME = datetime.datetime.now()

fig, axes = plt.subplots(1, 2)

colours = ['r', 'g', 'b', 'y', 'cyan', 'w', 'm']

axes[0].imshow(np.arcsinh(series[20]), cmap='viridis', interpolation='none', origin='lower')
axes[0].axis('off')
axes[0].set_title("Source Image")
for (xp, yp), c in zip(pos_1[:len(colours)], colours):
    circ = plt.Circle((xp, yp), 7, fill=False, edgecolor=c, linewidth=1)
    axes[0].add_patch(circ)

axes[1].imshow(np.arcsinh(series[0]), cmap='viridis', interpolation='none', origin='lower')
axes[1].axis('off')
axes[1].set_title("Target Image")
for (xp, yp), c in zip(pos_0[:len(colours)], colours):
    circ = plt.Circle((xp, yp), 7 * p.scale, fill=False, edgecolor=c, linewidth=1)
    axes[1].add_patch(circ)

plt.show()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')

# nickelodeon rhapsody
# crawling through africa 
# shrek parade