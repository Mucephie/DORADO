import warnings
warnings.filterwarnings('ignore')

from astropy.nddata.ccddata import CCDData
from ceres import Ceres, Stack
import numpy as np
from astropy import config as _config
import astropy.units as un
from ccdproc import ImageFileCollection
import ccdproc

import os
import datetime
from pathlib import Path




'''
Clippy is the handler of the Dorado system,
'''

# __all__ = []


# os.environ['HOME'] is also the home directory if needed
# os.uname() for returning operating system information

# os.open os.read os.write
# os.get_terminal_size (gets columns, lines)

# os.chdir for choosing current directory path
# os.getcwd() get the current working directory
# os.listdir(path='.') list directory contents by name

# os.mkdir(path) make a directory
# os.mkdirs(name, exist_ok = False) make directory and intermediate directories

# os.remove(path) remove/delete a file path (not for directories)
# os.rmdir(path) remove/delete a directory path (if dir is empty)
# os.removedirs(name) recursively removes directories until it finds a non empty directory

# os.rename(src, dst) rename a file or directory from src to dst
# os.scandir(path='.') better than listdir, you can specify a file descriptor
# with os.scandir(path) as it:
#     for entry in it:
#         if not entry.name.startswith('.') and entry.is_file():
#             print(entry.name)

# os.startfile(path[, operation]) acts like double clicking file in explorer to open with default
# os.system(command) execute a system command string in subshell

# getpass.getuser() similar to os.getlogin() gets username or login name. 




class clippy:

    def __init__(self):
        # open and use logger
        # make function to create data class from hardware, processed, or raw data folder
        # needs function to import data into raw
        # fix get_night()
        # add new calibration frames to folders
        # find most recetn calibration frame if none
        # -> figure out if none
        # clear cache function
        self.stardir = os.getcwd()
        self.rootname = 'dorado'
        self.config_dir = _config.get_config_dir(self.rootname)
        self.dordir = Path(self.config_dir).parent
        self.init_dir()
        self.unit = un.adu

    def init_dir(self):
        self.enter_dordir()
        os.makedirs('../data/wrk', exist_ok = True)
        os.makedirs('../data/flats', exist_ok = True)
        os.makedirs('../data/bias', exist_ok = True)
        os.makedirs('../data/darks', exist_ok = True)
        os.makedirs('../data/raw', exist_ok = True)
        os.makedirs('../data/graphical', exist_ok = True)
        os.makedirs('../data/projects', exist_ok = True)
        os.makedirs('../targets', exist_ok = True)
        os.makedirs('../logs', exist_ok = True)
        os.makedirs('../cache', exist_ok = True)

    def newdat(self):
        # find data that hasn't been processed yet
        print('searching for unprocessed data...')

    def get_night(self):
        """
        get_night obtains a timestring for the most recent(previous) night based on local/hardware
        time provided by datetime. The format follows yyyy-mm-(dd-1)+dd where (dd-1) is last nights 
        day of the month. This currently does not support dates at which the start of the observing
        night was the last day of the month.

        Parameters
        ----------
        None

        Returns
        -------
        night: str
                Timestring for the most recent night.
        """
        # currently does not support first/last of the month
        # option for last night or tonight
        date = datetime.date.today()
        year = date.year
        month = date.month
        date2 = date.day
        date1 = date2 - 1
        night = str(year) + '-' + str(month) + '-' + str(date1) + '+' + str(date2)

        return night

    def enter_dordir(self):
        os.chdir(self.dordir)

    def exit_dordir(self):
        os.chdir(self.stardir)

    def diread(self, date):
        path = self.dordir / 'data' / 'raw' / date
        contents = os.scandir(path = path)
        files = []
        directories = []
        for entry in contents:
            if not entry.name.startswith('.'):
                if entry.is_file():
                    files.append(entry)
                if entry.is_dir():
                    directories.append(entry)
        return files, directories
    
    def dirscan(self, date):

        files, directories = self.diread(date)

        biasstr = ['Bias', 'Bias', 'bias']
        flatsstr = ['FLAT', 'FlatField', 'flat', 'Flat', 'Flats', 'flats', 'FLATS', 'FlatFields']
        lightsstr = ['lights', 'Lights', 'LIGHTS']

        if len(directories) == 0:
            if len(files) == 0:
                raise Exception('No viable data found')

            else:
                print('Single directory level organization format detected.')
                print('Reading files.')
                # compile these into master frame and pass them to clippy
                path = self.dordir / 'data' / 'raw' / date
                # if np.any(for i in biasstr if i in s.name)
                biasl = [s for s in files if s.name in biasstr]
                bias = ImageFileCollection(location = path, filenames = biasl)
                flatsl = [s for s in files if s.name in flatsstr]
                flats = ImageFileCollection(location = path, filenames = flatsl)
                # strip into ceres (check if multi-filter)
                lightsl = [s for s in files if (s.name not in flatsstr) and (s.name not in biasstr)]
                lights = ImageFileCollection(location = path, filenames = lightsl)
                
                return bias, flats, lights


        elif len(files) == 0:
            print('Multi directory level organization format detected.')
            # print('Reading directories.')
            # read into this and compile into master bias, pass to clippy
            biasdir = [s for s in directories if s.name in biasstr]
            # check if multifilter, compile, pass to clippy
            flatsdir = [s for s in directories if s.name in flatsstr]
            # check if multifilter (or subdirectories) and pass to ceres
            lightsdir = [s for s in directories if (s.name not in flatsstr) and (s.name not in biasstr)]
            biasl, _ = self.diread(biasdir[0])
            bias = []
            for i in biasl:
                hdu = CCDData.read(i.path, unit = self.unit)
                bias.append(hdu)


            for ldir in lightsdir:
                files, directories = self.diread(ldir)
                if len(directories) == 0:
                    if len(files) == 0:
                        raise Exception('No viable light data found')
                    else:
                        print('Single directory lights organization format detected.')
                        # print('Reading files.')
                        lightsarr = []
                        for i in files:
                            hdu = CCDData.read(i.path, unit = self.unit)
                            lightsarr.append(hdu)
                        filter = ImageFileCollection(ldir).values('filter', unique = True)
                        lights = [filter, lightsarr]

                        
                elif len(files) == 0:
                    print('Multi directory lights organization format detected.')
                    # print('Reading directories.')
                    lights = []



            for fdir in flatsdir:
                files, directories = self.diread(fdir)
                if len(directories) == 0:
                    if len(files) == 0:
                        raise Exception('No viable light data found')
                    else:
                        print('Single directory flats organization format detected.')
                        # print('Reading files.')
                        flats = []
                        for i in files:
                            hdu = CCDData.read(i.path, unit = self.unit)
                            flats.append(hdu)

                elif len(files) == 0:
                    print('Multi directory flats organization format detected.')
                    # print('Reading directories.')
                    flats = []

           
            return bias, flats, lights

    def mkFlat(self, flats):
            """
            mkFlat takes  a list of flats to construct a calibrated flatfield image.
            
            Parameters
            ----------
            flats: array[CCDdata]
                    array of raw flatfields. ------------> this needs to be corrected for the new image storage format

            Returns
            -------
            flat: CCDdata
                    The combined calibrated flatfield image.
            """
            c = ccdproc.Combiner(flats)
            c.sigma_clipping()
            flat = c.median_combine()
            # , method = 'average',
            #                     sigma_clip = True, sigma_clip_low_thresh = 5, sigma_clip_high_thresh = 5,
            #                     sigma_clip_func = np.ma.median, sigma_clip_dev_func = mad_std, unit = self.unit)
            flat.header['stacked'] = True
            flat.header['numsubs'] = len(flats)

            return flat

    def mkBias(self, biasIFC):
            """
            mkBias takes a list of bias images to construct 
            a combined bias image. ------------> this needs to be corrected for the new image storage format
            Parameters
            ----------
            biasIFC: array[CCDdata]
                    array of raw bias images.

            Returns
            -------
            bias: CCDdata
                    The combined bias image.
            """
            # Allow specification of median or mean

            bias = ccdproc.combine(biasIFC, method = 'average', unit = self.unit)
            bias.meta['stacked'] = True
            bias.header['numsubs'] = len(biasIFC)

            return bias

    def mkceres(self, date):
            biasIFC, flats, lights = clip.dirscan(date)
            flat = self.mkFlat(flats)
            bias = self.mkBias(biasIFC)
            # save these frames
            time = [] # get_time() function

            filter = lights[0][0]
            lights = lights[1]
            print(filter)

            cere = Ceres()
            cere.bias = bias
            cere.flats = flat
            cere.add_stack(Stack(lights, filter = filter, calibrated = False), filter)

            return cere



clip = clippy()

# 2020-11-04+05
cere = clip.mkceres('2020-11-04+05')


# save the calibration frames to dorado's directory (check if that one already exists)
# we can compare these over time to look for its deviance

# is this calibrated data? aligned?
# extract the timestamps from each image

# we need to pass the calibration files(if none grab the most recent out of dorado by closest mjd)
# if there is no calibration files in dorado send a warning and set calibration files to none
