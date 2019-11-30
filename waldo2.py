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
from astropy.stats import mad_std
import astroalign
from astropy.time import Time

## photometry imports
from astropy import units as u
from photutils.psf import extract_stars
from photutils import EPSFBuilder

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

colours = ['r', 'g', 'b', 'y', 'cyan', 'w', 'm']

stardir = directory + night + '/wrk/stars/'
os.chdir(stardir)

## changes
# draco.plt_fits(series[3], 'viridis')
star3 = []
# # bad images
exclusions = [] # [3, 4, 5, 6, 91]
print('\nStarting star finding...\n')

# for q in range(80, len(series)):
#         if (q in exclusions):
#                 continue
#         data = series[q]
#         starT, star4 = w2f.dypeg_waldo(data)
#         if (starT == 0):
#                 print('\nToo little stars found, exclude ', q)
#         else:
#                 starT.write(str(q) + '-' + target + 'stars.fits')
#                 star3.append(star4)
#         if (q%5 == 0):
#                 print(q)
#                 print('\n', (q/5)/((len(series)-len(series)%5)/5), '% completed')

# data = series[30]
# x, y, r, opt_img = dracoOP2.starSeeker2(data)


# print('\n')
# size = 100
# hsize = (size - 1) / 2

# mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)))

# stars_tbl = Table()
# stars_tbl['x'] = x[mask]
# stars_tbl['y'] = y[mask]
# nddata = NDData(data=data)
# starsq = extract_stars(nddata, stars_tbl, size=50)

# l = []
# for u in range(len(x[mask])):
#         l.append(np.sum(starsq[u]))

# stars = np.zeros((len(x[mask]), 4))
# stars[:, 0] = x[mask]
# stars[:, 1] = y[mask]
# stars[:, 2] = r[mask]
# stars[:, 3] = l


# stars = stars[stars[:,3].argsort()]
# stars = np.flip(stars, axis = 0)

# star5 = stars[(0,1),:]

# dracoOP2.plt_stars(data, star5[:, 0], star5[:, 1], star5[:, 2])







print('\nReading...')
curdir = os.getcwd()
dirlist = os.listdir(curdir)
for h in range(len(dirlist)):
        star3.append(Table.read(dirlist[h]))

g = []
for d in range(len(series)):
        if (d in exclusions):
                 continue
        g.append(d)

print('\nPlotting...')

for w in range(0, len(star3), 5):
        star5 = star3[w]
        data = series[g[w]]
        plt.style.use(astropy_mpl_style)
        plt.figure()
        up, down = dracoOP2.plt_eye(data)
        plt.imshow(data, cmap='viridis', vmin=down, vmax=up)
        cbar = plt.colorbar()
        plt.grid(False)
        colours = ['r', 'g', 'b', 'y', 'cyan', 'w', 'm']
        plt.scatter(star5['x'].data, star5['y'].data, s = star5['r'].data ** 2, edgecolors = colours, alpha = 0.5) 

        labels = ['DY Peg', 'GSC 0712 00542', '2MASS 23084645+1715182']

        for t in range(len(star5['x'].data)):
                plt.annotate('{}'.format(labels[t]),
                        xy=(star5['x'][t], star5['y'][t]), xycoords='data',
                        xytext=(-150, 130), textcoords='offset points', size='16',
                        arrowprops=dict(arrowstyle="->"))
        plt.show()

timetaken = []
starA = []
starB = []
diff = []
for w in range(len(star3)):
        star5 = star3[w]
        data = series[g[w]]
        datet = data.header['DATE-OBS'].to_datetime()
        timetaken.append(str(datet.hour) + ':' + str(datet.minute))
        stars_tbl = Table()
        stars_tbl['x'] = star5['x']
        stars_tbl['y'] = star5['y']
        nddata = NDData(data=data)
        stars = extract_stars(nddata, stars_tbl, size=50)
        # print(w)
        if (w in [3, 4, 8, 9, 27]):
                diff.append(0)
        else:
                print(w, '\n')
                epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, progress_bar=True)  
                epsf, fitted_stars = epsf_builder(stars)
                starA.append(np.sum(fitted_stars[0]))
                starB.append(np.sum(fitted_stars[1]))
                diff.append(np.sum(fitted_stars[0]) - np.sum(fitted_stars[1]))

t = Time(timetaken, format='fits')

plt.figure()

plt.plot(t.value, diff)

plt.show()


## end of changes

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')

# >>> from astropy.time import Time
# >>> from astropy.coordinates import SkyCoord
# >>> tm = Time(['2000:002', '2002:345'])
# >>> sc = SkyCoord([10, 20], [-45, +40], unit='deg')
# >>> t = Table([tm, sc], names=['time', 'skycoord'])
# >>> t
