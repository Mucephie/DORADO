import warnings
warnings.filterwarnings('ignore')

from ..ceres import Ceres
from ..stack import Stack
from ..core.readerClass import *
from ..dorphot.dorphotClass import *
from ..target.targetClass import *

import ccdprocx

from astropy.nddata.ccddata import CCDData
import numpy as np
from astropy import config as _config
from astropy.utils.misc import isiterable
import astropy.units as un
from astropy.time import Time


from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions  import TimeoutError
from astropy.wcs import WCS

import os
import datetime
from pathlib import Path


ast = AstrometryNet()

'''
Filer is the handler of the Dorado system,
'''

__all__ = ['Dorado']

 
def get_builtins():
    """Due to the way Python works, ``__builtins__`` can strangely be either a module or a dictionary,
    depending on whether the file is executed directly or as an import. I couldn’t care less about this
    detail, so here is a method that simply returns the namespace as a dictionary."""
    return getattr( __builtins__, '__dict__', __builtins__ )


class Dorado_core:
    '''
    
    '''
    def __init__(self):
        # open and use logger
        # make function to create data class from hardware, processed, or raw data folder
        # needs function to import data into raw
        # fix get_night()
        # add new calibration frames to folders
        # find most recent calibration frame if none
        # -> figure out if none
        # clear cache function
        self.stardir = os.getcwd()
        self.rootname = 'dorado'
        self.config_dir = _config.get_config_dir(self.rootname)
        self.dordir = Path(self.config_dir).parent
        self.init_dir()
        self.unit = un.adu

        ## class storage zone -----------------------
        # reader
        self.reader = aico_reader() # ussually we'd want to grab this from a configuration file, should name be aicoReader?
        # dorphot
        self.dorphot = aicoPhot()
        # ceres & stacks
        self.ceres_keys = {}
        self.ceres = [] # maybe call this series?
        # self.filters[stack.filter] = len(self.data)
        # self.data.append(stack)

        # targets
        self.target_keys = {}
        self.targets = [] # TODO choose target class
        # all early ts instances will belong to a target;
        # timeseries
        self.timeseries = [] # exclude this for the time being as 

        # site information (hardcoded atm)
        self.UTCoffset = -5
        
    def init_dir(self, tess = False):
        """
        
        Parameters
        ----------
        tess: boolean
            Whether to initialize a TESS data directory in './data/'
        Returns
        -------
        
        """
        self.enter_dordir()
        os.makedirs('./data/wrk', exist_ok = True)
        os.makedirs('./data/flats', exist_ok = True)
        os.makedirs('./data/bias', exist_ok = True)
        os.makedirs('./data/darks', exist_ok = True)
        os.makedirs('./data/raw', exist_ok = True)
        os.makedirs('./data/graphical', exist_ok = True)
        os.makedirs('./data/projects', exist_ok = True)
        os.makedirs('./data/targets', exist_ok = True)
        os.makedirs('./logs', exist_ok = True)
        os.makedirs('./cache', exist_ok = True)
        os.makedirs('./cache/astrometryNet', exist_ok = True)
        if tess:
            os.makedirs('./data/tess', exist_ok = True)
        self.exit_dordir()
        
    def newdat(self):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
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
        # UTC or local time?
        date = datetime.date.today()
        year = date.year
        month = date.month
        day = date.day
        # date1 = date2 - 1
        # night = str(year) + '-' + str(month) + '-' + str(date1) + '+' + str(date2)


        if day < 10:
            daystr = '0' + str(day)
        else:
            daystr = str(day)
        if day < 9:
            day2str = '0' + str(day + 1)
        else:
            day2str = str(day + 1)
        if month < 10:
            monthstr = '0' + str(month)
        else:
            monthstr = str(month)


        night = str(year) + '-' + monthstr + '-' + daystr + '+' + day2str




        return night
        
    def enter_dordir(self):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        os.chdir(self.dordir)
         
    def exit_dordir(self):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        os.chdir(self.stardir)
        
    def diread(self, dirarray):
        """
        
        Parameters
        ----------
        dirarray: str array
            array containing the path to the desired directory.
        Returns
        -------
        
        """
        if isiterable(dirarray):
            path = self.dordir
            for dir in dirarray:
                path = path / dir
        else:
            path = dirarray
        # path = self.dordir / 'data' / 'raw' / date
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
    
    def mkceres(self, date, name = None, sub = 'raw', target = None, calibrated = False, aligned = False):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        if target == None:
            ctemp = self.reader.mkceres(date, sub = sub, target = target, calibrated = calibrated, aligned = aligned)
        else:
            ctemp = self.reader.mkceres(date, sub = sub, target = self.targets[self.target_keys[target]], calibrated = calibrated, aligned = aligned)

        if name == None:
            name = ctemp.datestr
            # TODO handle if no date found either
            print('No series nickname given, defaulting to date string: ', name)
        else:
            print('Call series as: ', name)
        self.ceres_keys[name] = len(self.ceres)
        self.ceres.append(ctemp)
    
    def mktrgt(self, name):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        # add dictionary
        # , coordinates = None
        self.target_keys[name] = len(self.targets)
        self.targets.append(Target(name))
        
    def force16(self, hdu):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        print('This function is not implemented yet. See mkBias() for example functionality.')
        
    def mkcacheObj(self, object, subcache = False):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        if subcache:
            cachedir = self.dordir / 'cache' / subcache
            dirarray = ['cache', subcache]
        else:
            cachedir = self.dordir / 'cache' 
            dirarray = ['cache']
        # if isiterable(object):
        #     # print('This function does not currently support iterable objects.')
        #     return warnings.WarningMessage('This function does not currently support iterable objects.')
        # else:
        #     files, _ = self.diread(dirarray)
        #     fname = 'cache_object_' + str(len(files) + 1) + '.fits'
        #     object.write(fname)
        #     return(fname, cachedir)

        # check if iterable
        files, _ = self.diread(dirarray)
        fname = 'cache_object_' + str(len(files) + 1) + '.fits'
        
        object.write(cachedir / fname, overwrite = True)
        return(fname, cachedir)
        
    def delcacheObj(self, fname, subcache = False):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        # TODO extend to clearing cache
        if subcache:
            cachedir = self.dordir / 'cache' / subcache
        else:
            cachedir = self.dordir / 'cache'

        os.remove(cachedir / fname)
        
    def plate_solve(self, dirarray, data = None, writearray = False):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        path = self.dordir
        for dir in dirarray:
            path = path / dir

        if data == None:
            data = CCDData.read(path) #, unit = Dorado.unit) ## NOTE edited

        trying = True
        submission_id = None
        num = 0

        while trying:
                try:
                    if not submission_id:
                        wcs_header = ast.solve_from_image(path, force_image_upload=True, submission_id=submission_id, solve_timeout=300)
                    else:
                        print('Monitoring: try #', num)
                        wcs_header = ast.monitor_submission(submission_id, solve_timeout=300)
                except TimeoutError as e:
                    print(TimeoutError)
                    num = num + 1
                    print('Timed out: try #', num)
                    submission_id = e.args[1]

                if wcs_header != None:
                    # got a result, so terminate while loop
                    trying = False
        if wcs_header:
            # Code to execute when solve succeeds
            print('Solve succeeded! :)')
            wcs_hdu = data
            wcs_hdu.header = wcs_header
            if writearray:
                path = self.dordir
                for dir in writearray:
                    path = path / dir
                wcs_hdu.write(path, overwrite = True)

            return wcs_hdu, wcs_header
        else:
            # Code to execute when solve fails
            print('Solve failed! :(')
            return 
         
    def getDateString(self, cr):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """

        ## TODO :: look into UTC wrecking stuff
        # if the hour is less than the utc offset of the site then the utc date is one ahead of local time
        epoch = self.ceres[self.ceres_keys[cr]].date.ymdhms

        # epoch = date.ymdhms

        if (epoch['hour'] + self.UTCoffset) < 0:
            day = str(epoch['day'] - 1)
            day2 = str(epoch['day'])
            month = str(epoch['month'])

            if (epoch['day'] - 1) < 10:
                day = '0' + str(epoch['day'] - 1)
            if epoch['day'] < 10:
                day2 = '0' + str(epoch['day'])

            if epoch['month'] < 10:
                month = '0' + str(epoch['month'])

        else: 
            day = str(epoch['day'])
            day2 = str(epoch['day'] + 1)
            month = str(epoch['month'])

            if epoch['day'] < 10:
                day = '0' + str(epoch['day'])
                if epoch['day'] < 9:
                    day2 = '0' + str(epoch['day'] + 1)

            if epoch['month'] < 10:
                month = '0' + str(epoch['month'])

        datestr = str(epoch['year']) + '-' + month + '-' + day + '+' + day2
        self.ceres[self.ceres_keys[cr]].datestr = datestr
        
    def mkwrk(self, cr):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        # TODO this should be a wrapper for reader
        if self.ceres[self.ceres_keys[cr]].datestr == None:
            self.getDateString(cr) 
        datestr = self.ceres[self.ceres_keys[cr]].datestr
        wrkdir = self.dordir / 'data' / 'wrk'
        os.makedirs(wrkdir / datestr, exist_ok = True)
        os.makedirs(wrkdir / datestr / 'aligned', exist_ok = True)
        os.makedirs(wrkdir / datestr / 'calibrated', exist_ok = True)
        os.makedirs(wrkdir / datestr / 'uncalibrated', exist_ok = True)
        os.makedirs(wrkdir / datestr / 'WCS', exist_ok = True)
        # figures, targets, log, omitted images, observation metadata
        
    def savewrk(self, cr, filters = None):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        self.reader.savewrk(cr, filters)
        
    def saveWCS(self, cr, filters = None):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        wrkdir = self.dordir / 'data' / 'wrk'
        if self.ceres[self.ceres_keys[cr]].datestr == None:
            self.getDateString(cr)
        datestr = self.ceres[self.ceres_keys[cr]].datestr
        self.mkwrk(cr)
        if filters == None:
            filters = self.ceres[self.ceres_keys[cr]].filters.keys()
        
        for filter in filters:
            stack = self.ceres[self.ceres_keys[cr]].data[self.ceres[self.ceres_keys[cr]].filters[filter]]
            if stack.wcs == None:
                self.dorphot.getWCS(filter, self) # TODO this is now in dorphot, should it be moved?
            fname = str(filter) + '-solved.fits'
            solved = self.ceres[self.ceres_keys[cr]].data[self.ceres[self.ceres_keys[cr]].filters[filter]].solved
            solved.write(wrkdir / datestr / 'WCS' / fname, overwrite = True)
        # not done might be redundant atm also shouldnt this belong to reader?
    def saveBase(self, cr, filters = None):
        """
        
        Parameters
        ----------

        Returns
        -------
        
        """
        wrkdir = self.dordir / 'data' / 'wrk'
        if self.ceres[self.ceres_keys[cr]].datestr == None:
            self.getDateString(cr)
        datestr = self.ceres[self.ceres_keys[cr]].datestr

        print('Saved to data/wrk/', datestr)
        self.mkwrk(cr)

        if filters == None:
            filters = self.ceres[self.ceres_keys[cr]].filters.keys()
        
        for filter in filters:
            imname = filter + '_base.fits'
            fname = wrkdir / datestr / imname
            base = self.ceres[self.ceres_keys[cr]].data[self.ceres[self.ceres_keys[cr]].filters[filter]].base
            base.write(fname, overwrite = True)
        # not done might be redundant atm also shouldnt this belong to reader?
        

    # merge header

Dorado = Dorado_core()

# this gone pollute the namespace but ay lmao
G = get_builtins()
G[ 'G' ] = G
G['Dorado'] = Dorado



# we need to pass the calibration files(if none grab the most recent out of 
# dorado by closest mjd)
# if there is no calibration files in dorado send a warning and set 
# calibration files to none

# make a function to force uint16




# if aligned:
            #     dirarray = ['data', sub, date, 'aligned']
            # elif calibrated:
            #     dirarray = ['data', sub, date, 'calibrated']
            # else:
            #     dirarray = ['data', sub, date]
            # biasIFC, flats, lights = self.dirscan(dirarray)
            # print(len(flats), ' flats found.')
            # print(len(biasIFC), ' bias frames found.')
            # print(len(lights), ' lights found')
            
            

            # # save these frames
            # ## TODO :: get missing calibraation files

            # ## TODO :: look into UTC wrecking stuff
            # if len(biasIFC) == 0:
            #     cere = Ceres(time = Time(lights[0].header['DATE-OBS'], format='fits'))
            #     self.getDateString(cere)
            # else:
            #     bias = self.mkBias(biasIFC)
            #     cere = Ceres(bias = bias, time = Time(lights[0].header['DATE-OBS'], format='fits'))
            #     self.getDateString(cere)

            # # for f in filter_data:
            # # cere.flats[flat.header['filter']] = flat

            # ## TODO :: multifilter fun
            # if len(flats) == 0:
            #      cere.add_stack(Stack(lights, calibrated = calibrated, aligned = aligned, target = target))
            # else:
            #     flat = self.mkFlat(flats)
            #     cere.add_stack(Stack(lights, flat = flat, calibrated = calibrated, aligned = aligned, target = target))
            

            # return cere
