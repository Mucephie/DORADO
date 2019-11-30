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
import skimage.exposure as skie
from scipy.ndimage import median_filter
from skimage.feature import blob_dog, blob_log
import astroalign
from scipy.ndimage import zoom

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




## functions
def mkwrkdir(datadir, night, expath = ''):
        stardir = os.getcwd()

        path = datadir + night + '/' + expath # + '/'
        os.chdir(path)
        os.mkdir(path + '/wrk')
        path = path + '/wrk'
        os.mkdir(path + '/bias')
        os.mkdir(path + '/flats')
        os.mkdir(path + '/lights')
        os.mkdir(path + '/calibrated')
        os.mkdir(path + '/stars')
        os.mkdir(path + '/results')

        os.chdir(stardir)

def get_night():
      # currently does not support first/last of the month
      date = datetime.date.today()
      year = date.year
      month = date.month
      date2 = date.day
      date1 = date2 - 1
      night = str(year) + '-' + str(month) + '-' + str(date1) + '+' + str(date2)
      return night

def checkdir(directory, night, expath = ''):
        path = directory + '/' + night + '/' + expath
        print(path)
        dirlist = os.listdir(path)
        bias = [s for s in dirlist if 'BIAS' in s]
        if len(bias)==0:
                bias = [s for s in dirlist if 'Bias' in s]
        if len(bias)==0:
                bias = [s for s in dirlist if 'Bias' in s]
        flats = [s for s in dirlist if 'FLAT' in s]
        if len(flats)==0:
                flats = [s for s in dirlist if 'FlatField' in s]
        lights = [s for s in dirlist if (np.invert('FLAT' in s)) and (np.invert('Flat' in s)) and (np.invert('BIAS' in s)) and (np.invert('Bias' in s)) and ('.FIT' in s) ]
        if len(lights)==0:
                lights = [s for s in dirlist if (np.invert('FLAT' in s)) and (np.invert('Flat' in s)) and (np.invert('BIAS' in s)) and (np.invert('Bias' in s)) and ('.fit' in s) ]
        if len(lights)==0:
                lights = [s for s in dirlist if (np.invert('FLAT' in s)) and (np.invert('Flat' in s)) and (np.invert('BIAS' in s)) and (np.invert('Bias' in s)) and ('.fits' in s) ]
        if len(lights)==0:
                lights = [s for s in dirlist if (np.invert('FLAT' in s)) and (np.invert('Flat' in s)) and (np.invert('BIAS' in s)) and (np.invert('Bias' in s)) and ('.FITS' in s) ]
        print('\ndirlist: ', len(dirlist))
        print('\nbias\': ', len(bias))
        print('\nflats: ', len(flats))
        print('\nlights: ', len(lights))
        print('\ntotal: ', len(lights) + len(flats) + len(bias))

        return bias, flats, lights

def get_series(directory, night, imlist, expath = '', unit = u.adu):
        ## rewrite to match get series in reduce_series
        stardir = os.getcwd()
        # print(stardir)
        path = directory + night + '/' + expath
        # print(path)
        os.chdir(path)
        series = []
        for i in range(len(imlist)):
                os.chdir(path)
                hdu = CCDData.read(imlist[i], unit=unit)
                series.append(hdu)
        # series = ccdproc.ImageFileCollection(filenames=imlist)
        os.chdir(stardir)


        return series

def write_series(directory, night, series, vstar):
        stardir = os.getcwd()
        caldir = directory + night + '/wrk/calibrated'
        os.chdir(caldir)
        for i in range(len(series)):
                series[i].write(vstar + str(i) + '.calibrated.fit')


        os.chdir(stardir)

def mastBias(directory, night, bias):
        stardir = os.getcwd()

        path = directory + '/' + night + '/'
        os.chdir(path)

        master_bias = ccdproc.combine(bias, method='average', unit=u.adu)
        master_bias.meta['stacked'] = True


        caldir = directory + night + '/wrk/bias/'
        os.chdir(caldir)
        master_bias.write('master_bias.fit')

        os.chdir(stardir)

        return master_bias

def mastFlat(directory, night, flats, bias):
        stardir = os.getcwd()

        path = directory + '/' + night + '/'
        os.chdir(path)
        # for hdu, i in flats.hdus():
        #         hdu = ccdproc.subtract_bias(hdu, bias)
        #         hdu.header['bias corrected'] = True
        #         hdu.header['calibrated'] = True
                

        master_flat = ccdproc.combine(flats, method='average',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, unit=u.adu)
        master_flat.header['stacked'] = True


        caldir = directory + night + '/wrk/flats/'
        os.chdir(caldir)
        master_flat.write('master_flat.fit')

        os.chdir(stardir)

        return master_flat

def reduce_series(directory, night, imlist, flat, bias, vstar, expath = ''):
        stardir = os.getcwd()
        imdir = directory + '/' + night + '/' + expath
        caldir = directory + night + '/wrk/calibrated/'
        os.chdir(caldir)
        series = []
        for i in range(len(imlist)):
                os.chdir(imdir)
                hdu = CCDData.read(imlist[i], unit=u.adu)
                hdu = ccdproc.ccd_process(hdu, master_bias = bias, master_flat = flat)
                hdu = ccdproc.cosmicray_lacosmic(hdu, sigclip=5)
                hdu.header['bias corrected'] = True
                hdu.header['flat corrected'] = True
                hdu.header['cosmicray corrected'] = True
                hdu.header['calibrated'] = True
                os.chdir(caldir)
                hdu.write(vstar + expath + '-' + str(i) + '-calibrated.fit')
                series.append(hdu)
        # print(len(series))
        os.chdir(stardir)

        return series

def reduce_seriesR(directory, night, imlist, flat, bias, vstar, resize,expath = ''):
        stardir = os.getcwd()
        imdir = directory + '/' + night + '/' + expath
        caldir = directory + night + '/wrk/calibrated/'
        os.chdir(caldir)
        series = []
        for i in range(len(imlist)):
                os.chdir(imdir)
                hdu = CCDData.read(imlist[i], unit=u.adu)
                hdu = ccdproc.ccd_process(hdu, master_bias = bias, master_flat = flat)
                hdu = ccdproc.cosmicray_lacosmic(hdu, sigclip=5)
                hdu.header['bias corrected'] = True
                hdu.header['flat corrected'] = True
                hdu.header['cosmicray corrected'] = True
                hdu.header['calibrated'] = True
                hdu.header['Resized'] = True
                hdu.data = zoom(hdu.data, (resize, resize), order=0)
                os.chdir(caldir)
                hdu.write(vstar + expath + '-' + str(i) + '-calibrated.fit')
                series.append(hdu)
        # print(len(series))
        os.chdir(stardir)

        return series

def plt_stars(data, x, y, r):
        plt.style.use(astropy_mpl_style)
        plt.figure()
        # plt.imshow(data, cmap='viridis', vmin=0)
        up, down = plt_eye(data)
        plt.imshow(data, cmap='viridis', vmin=down, vmax=up)
        cbar = plt.colorbar()
        plt.grid(False)
        
        # for i in range(len(x)):
        #         circlei=plt.Circle((x[i],y[i]), r[i], edgecolor=, alpha = 0.75, linewidth = 1)
        #         plt.gcf().gca().add_artist(circlei)
        colours = ['r', 'g', 'b', 'y', 'cyan', 'w', 'm']
        plt.scatter(x, y, s = r ** 2, edgecolors = colours, alpha = 0.5) 
        for p in range(len(x)):
                plt.text(x[p], y[p], str(p))

        plt.show()


def starSeeker2(data):
        mf = median_filter(data, size= 15)
        datamf = data - mf
        limg = np.arcsinh(datamf) #datamf
        limg = limg / limg.max()
        low = np.percentile(limg, 0.2)
        high = np.percentile(limg, 99.1)
        opt_img  = skie.exposure.rescale_intensity(limg, in_range=(low,high))

        stars =  blob_log(opt_img, max_sigma=100, min_sigma = 5, num_sigma=10, threshold=.2)
        # stars =  blob_dog(opt_img, max_sigma=30, threshold=.2)

        # Compute radii in the 3rd column.
        stars[:, 2] = stars[:, 2] * np.sqrt(2)

        y2, x2, r = stars[:, 0], stars[:, 1], stars[:, 2]


        print('\nStars found: ', len(r))

        return x2, y2, r, opt_img

def plt_eye(data):
    mean = np.mean(data)
    std = np.std(data)
    mean_up = mean + 2 * std
    mean_down = mean - 2 * std

    return mean_up, mean_down

def starange(x, y, r, field, fieldinfo, spread = 5):
        star = -1 * np.zeros(x.shape)
        isStar = []
        for k in range(len(x)): # Stars found
                distances = []
                for j in range(len(x)): # other staars found
                        dist = np.sqrt((x[k]-x[j])**2 + (y[k]-y[j])**2)
                        distances.append(dist)
                for l in range(len(field['star1'].data)): # all field stars for identification
                        starmatch = -1
                        starQ = field['star' + str(l + 1)].data
                        for t in range(len(distances)): # over every distance
                                for i in range(len(starQ.data)): # star distances for star[l]
                                        starpos = starQ[i]
                                        # print(np.around(starpos - spread))
                                        if (((starpos - spread) < distances[t]) and (distances[t] < (starpos + spread))):
                                                starmatch = starmatch +1
                                                print('\nDiff ', np.around(np.abs(starpos - distances[t])))
                        if (starmatch > 2) :
                                star[k] = l
                                print('\nStar ', k, ' is found to be star ', l + 1, ' with a match # of ', starmatch)
                                isStar.append(k)

                if ((star[k]) < 0): 
                        print('\nStar ', k, ' was not identified!')

                if ((star[k] > 2)):
                        isStar.append(k)
                print('\nDone star ', k)
        xs = []
        ys = []
        for w in range(len(isStar)):
                xs.append(x[w])
                ys.append(x[w])

        fieldres = Table(fieldinfo['field#'])
        fieldres['star#'] = [isStar]
        fieldres['x'] = [xs]
        fieldres['y'] = [ys]
        fieldres['ra'] = [fieldinfo['ra']]
        fieldres['dec'] = [fieldinfo['dec']]
        fieldres['isV'] = [fieldinfo['isV']]
        fieldres['isC'] = [fieldinfo['isC']]

        return fieldres, isStar

def theMask(data, lx, hx, ly, hy):
        mask = np.zeros(data.shape, dtype=bool)
        mask[lx:hx, ly:hy] = True
        

        return mask

def get_vstar(name):
        if (name=='DYPEG'):
                print('Star not initialized yet.')
        elif (name=='XXCYG'):
                print('Star not initialized yet.')
        elif (name=='YZBOO'):
                print()
        elif (name=='BLCAM'):
                print('Star not initialized yet.')
        elif (name=='BELYN'):
                print('Star not initialized yet.')






# ## testing
# START_DATE_TIME = datetime.datetime.now()
# print('\nStarting time: ', START_DATE_TIME)

# home_dir = os.getcwd()
# print('\nHome dir: ', home_dir)
# target = input('Enter the target name (e.g DYPEG): ')
# usetdy = input('Reduce todays data (y/n)? ')
# if (usetdy == 'y'):
#         night = get_night()
# else:
#         night = input('Enter night code (YYYY-MM-DD+DD): ')
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
# directory = 'D:/draco_data/'
# mkwrkdir(directory, night)

# bias_list, flats_list, lights_list = checkdir(directory, night)

# biaseries = get_series(directory, night, bias_list)
# bias = mastBias(directory, night, bias_list)

# flats = get_series(directory, night, flats_list)
# flat = mastFlat(directory, night, flats_list, bias)

# series = get_series(directory, night, lights_list)
# series = reduce_series(directory, night, lights_list, flat, bias, target)

# END_DATE_TIME = datetime.datetime.now()

# print('\nEnding time: ', END_DATE_TIME)
# print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME), '\n')


