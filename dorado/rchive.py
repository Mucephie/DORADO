##  imports
import os
import datetime
from astropy.nddata import CCDData
from astropy.io import fits
import astropy.units as u
import numpy as np

__all__ = ['mkwrkdir', 'get_night', 'checkdir', 'get_series', 'write_series', 'get_im', 'file_string']

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
        path = directory + '/' + night + '/' + expath
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

def get_im(directory, imname):
    os.chdir(directory)
    file_string = imname
    reduc_data = fits.getdata(file_string)
    reduc_header = fits.getheader(file_string)

    return reduc_data, reduc_header

def file_string(directory, night, im, expath = ''):
        # this is a temporary function for ease of testing
        path = directory + '/' + night + '/' + expath + im
        return path

# def lightcurve_save
# def timing/logging
