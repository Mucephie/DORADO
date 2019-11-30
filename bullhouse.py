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
# from astropy.stats import mad_std
# import skimage.exposure as skie
# from scipy.ndimage import median_filter
# from skimage.feature import blob_dog, blob_log
import astroalign
from donuts import Donuts

## photometry imports
from astropy import units as u

## plotting imports
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.visualization import astropy_mpl_style
import matplotlib.colors as colors
from matplotlib import cm

## other imports
import datetime
import dracoOP2
import draco

## testing
START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)

# home_dir = os.getcwd()
# print('\nHome dir: ', home_dir)
# target = input('Enter the target name (e.g DYPEG): ')
# usetdy = input('Reduce todays data (y/n)? ')
# if (usetdy == 'y'):
#         night = dracoOP2.get_night()
# elif (usetdy == 'z'):
#         night = '2019-08-14+15'
# elif (usetdy == 'x'):
#         night = '2019-08-23+24'
# else:
#         night = input('Enter night code (YYYY-MM-DD+DD): ')

# # # data directory: D:/draco_data/
# # # test directory: C:/Users/mucep/Offline/Draco-test/
# # ## nights
# # # 2018-10-9+10 GALAXIES
# # # 2018-12-04+05 DYPEG $$
# # # 2018-12-18+19 46P & FSTAR
# # # 2019-01-04+05 BLCAM $$
# # # 2019-01-13+14 BLCAM $$
# # # 2019-08-14+15 YZBOO $$
# # # 2019-08-23+24 XXCYG $$
# directory = 'D:/draco_data/'
# datadir = directory + night + '/wrk/calibrated/'

# bias_list, flats_list, lights_list = dracoOP2.checkdir(directory, night, '/wrk/calibrated/')
# os.chdir(datadir)
# #series = dracoOP2.get_series(directory, night, lights_list, '/wrk/calibrated/', unit = None)

# ## bullhouse algorithm  
# refim = lights_list[0]

# d = Donuts(refimage=refim, image_ext=0, subtract_bkg=True, ntiles=64)

# # for each image, compute the x/y translation required
# # to align the images onto the reference image
# xp = []
# yp = []
# for image in lights_list:
#      shift_result = d.measure_shift(image)
#      x = shift_result.x
#      y = shift_result.y
#      xp.append(x.value)
#      yp.append(y.value)
#      # Also check out shift_result.sky_background
#      print(x, y)


# plt.figure()
# plt.plot(xp, yp, 'g-')
# plt.xlabel('x (px)')
# plt.ylabel('y (px)')
# plt.title('60cm Tracking pattern')

# plt.show()

# 1: 1005 748
# 50: 826 712
# 100: 655 652
# 131: 560 613


# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

xt = [1005, 918, 826, 763, 655,572, 560]
xt = np.array(xt) - xt[0]
yt = [748, 733, 712, 682, 652, 622, 613]
yt = np.array(yt) - yt[0]
zt = [1, 25, 50, 75, 100, 125, 131]

x2 = [1187, 1115, 1055, 1008, 963, 903, 853] 
x2 = np.array(x2) - x2[0]
y2 = [419, 407, 386, 363, 336, 311, 279] 
y2 = np.array(y2) - y2[0]
z2 = np.array([1, 25, 50, 75, 100, 125,150]) * 5/6

x3 = [1263, 1227, 1195, 1123, 1048, 987, 912] 
x3 = np.array(x3) - x3[0]
y3 = [669, 666, 651, 623, 590, 555, 515] 
y3 = np.array(y3) - y3[0]
z3 = np.array([1, 25, 50, 75, 100, 125,150]) * 2

x4 = [1367, 1281, 1206, 1117, 1059] 
x4 = np.array(x4) - x4[0]
y4 = [881, 845, 799, 748, 680] 
y4 = np.array(y4) - y4[0]
z4 = np.array([1, 25, 50, 75, 100]) * 2.5

x5 = [1496, 1457, 1371, 1277, 1220, 1148, 1102] 
x5 = np.array(x5) - x5[0]
y5 = [866, 858, 802, 759, 705, 652, 590] 
y5 = np.array(y5) - y5[0]
z5 = np.array([1, 25, 50, 75, 100, 125,150]) * 2

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', aspect='equal')
ax.plot(xt, yt, zt, '-')
ax.scatter(xt, yt, zt)
ax.plot(x2, y2, z2, '-')
ax.scatter(x2, y2, z2)
ax.plot(x3, y3, z3, '-')
ax.scatter(x3, y3, z3)
ax.plot(x4, y4, z4, '-')
ax.scatter(x4, y4, z4)
ax.plot(x5, y5, z5, '-')
ax.scatter(x5, y5, z5)

ax.set_xlabel('x (px)')
ax.set_ylabel('y (px)')
ax.set_zlabel('time (min)')
plt.title('60cm Tracking pattern')
plt.show()

## end of Bullhouse algorithm

## end of changes

END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')