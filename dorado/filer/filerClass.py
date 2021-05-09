import warnings
warnings.filterwarnings('ignore')

from ..ceres import Ceres
from ..stack import Stack


from astropy.nddata.ccddata import CCDData
import numpy as np
from astropy import config as _config
from astropy.utils.misc import isiterable
import astropy.units as un
from astropy.time import Time
import ccdproc

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

__all__ = ['Filer']

class Filer:

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
        self.exit_dordir()

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

    def diread(self, dirarray):
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
    
    def dirscan(self, dirarray):
        path = self.dordir
        for dir in dirarray:
            path = path / dir
        files, directories = self.diread(path)

        biasstr = ['Bias', 'Bias', 'bias', 'BIAS']
        flatsstr = ['FLAT', 'FlatField', 'flat', 'Flat', 'Flats', 'flats', 'FLATS', 'FlatFields']
        lightsstr = ['lights', 'Lights', 'LIGHTS']

        if len(directories) == 0:
            if len(files) == 0:
                raise Exception('No viable data found')

            else:
                print('Single directory level organization format detected.')
                # print('Reading files.')
                # compile these into master frame and pass them to Filer
                # path = self.dordir / 'data' / 'raw' / date
                biasl = []
                for strbias in biasstr:
                    for s in files:
                        if (str(strbias)) in str(s.name):
                            biasl.append(s)
                bias = []
                for i in range(len(biasl)):
                    hdu = CCDData.read(biasl[i].path) #, unit = self.unit)
                    bias.append(hdu)
                print('Bias searched.')
                flatsl = []
                for strflat in flatsstr:
                    for s in files:
                        if strflat in s.name:
                            flatsl.append(s)
                flats = []
                for i in flatsl:
                    hdu = CCDData.read(i.path) #, unit = self.unit)
                    flats.append(hdu)
                print('flats searched.')
                # strip into ceres (check if multi-filter)
                lightsl = []
                nonlight = []
                for b in biasl:
                    nonlight.append(b.name)
                for f in flatsl:
                    nonlight.append(f.name)

                for s in files:
                    if (str(s.name) not in nonlight):
                        lightsl.append(s)

                # lightsl = [s for s in files if (s.name not in biasl) and (s.name not in flatsl)]
                lights = []
                for i in lightsl:
                    hdu = CCDData.read(i.path) # , unit = self.unit)
                    lights.append(hdu)
                print('lights searched.')
                return bias, flats, lights


        elif len(files) == 0:
            print('Multi directory level organization format detected.')
            # print('Reading directories.')
            # read into this and compile into master bias, pass to Filer
            biasdir = [s for s in directories if s.name in biasstr]
            # check if multifilter, compile, pass to Filer
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
                        lights = []
                        for i in files:
                            hdu = CCDData.read(i.path, unit = self.unit)
                            lights.append(hdu)
                        # filter = ImageFileCollection(ldir).values('filter', unique = True)
                        # lights = [filter, lightsarr]

                        
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
            date = Time(bias.header['DATE-OBS'], format='fits').mjd
            bias.data = bias.data.astype('uint16') 
            fname = str(int(date)) + '_Bias.fits'
            biasdir = self.dordir / 'data' / 'bias' 
            contents = os.scandir(path = biasdir)
            save = True
            for entry in contents:
                if fname in entry.name:
                    save = False
            if save:
                print('Saving Bias for later use')
                bias.write(biasdir / fname)
            else:
                print('Bias for date already saved.')
            return bias

    def mkceres(self, date, sub = 'raw', target = None, calibrated = False, aligned = False):
            if aligned:
                dirarray = ['data', sub, date, 'aligned']
            elif calibrated:
                dirarray = ['data', sub, date, 'calibrated']
            else:
                dirarray = ['data', sub, date]
            biasIFC, flats, lights = self.dirscan(dirarray)
            print(len(flats), ' flats found.')
            print(len(biasIFC), ' bias frames found.')
            print(len(lights), ' lights found')
            
            

            # save these frames


            if len(biasIFC) == 0:
                cere = Ceres(time = Time(lights[0].header['DATE-OBS'], format='fits'))
            else:
                bias = self.mkBias(biasIFC)
                cere = Ceres(bias = bias, time = Time(lights[0].header['DATE-OBS'], format='fits'))

            # for f in filter_data:
            # cere.flats[flat.header['filter']] = flat

            if len(flats) == 0:
                 cere.add_stack(Stack(lights, calibrated = calibrated, aligned = aligned, target = target))
            else:
                flat = self.mkFlat(flats)
                cere.add_stack(Stack(lights, flat = flat, calibrated = calibrated, aligned = aligned, target = target))
            

            return cere

    def force16(self, hdu):
        print('This function is not implemented yet. See mkBias() for example functionality.')

    def mkcacheObj(self, object, subcache = False):
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
        if subcache:
            cachedir = self.dordir / 'cache' / subcache
        else:
            cachedir = self.dordir / 'cache'

        os.remove(cachedir / fname)
    
    def plate_solve(self, dirarray, data = None, writearray = False):
        path = self.dordir
        for dir in dirarray:
            path = path / dir

        if data == None:
            data = CCDData.read(path, unit = self.unit)

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
        day = str(cr.date.ymdhms['day'])
        day2 = str(cr.date.ymdhms['day'] + 1)
        month = str(cr.date.ymdhms['month'])
        if cr.date.ymdhms['day'] < 10:
            day = '0' + str(cr.date.ymdhms['day'])
            if cr.date.ymdhms['day'] < 9:
                day2 = '0' + str(cr.date.ymdhms['day'] + 1)
        if cr.date.ymdhms['month'] < 10:
            month = '0' + str(cr.date.ymdhms['month'])

        datestr = str(cr.date.ymdhms['year']) + '-' + month + '-' + day + '+' + day2
        cr.datestr = datestr

    def mkwrk(self, cr):
        if cr.datestr == None:
            self.getDateString(cr)
        datestr = cr.datestr
        wrkdir = self.dordir / 'data' / 'wrk'
        os.makedirs(wrkdir / datestr, exist_ok = True)
        os.makedirs(wrkdir / datestr / 'aligned', exist_ok = True)
        os.makedirs(wrkdir / datestr / 'calibrated', exist_ok = True)
        os.makedirs(wrkdir / datestr / 'uncalibrated', exist_ok = True)
        os.makedirs(wrkdir / datestr / 'WCS', exist_ok = True)
        # figures, targets, log, omitted images, observation metadata
        
    def savewrk(self, cr, filters = None):
        wrkdir = self.dordir / 'data' / 'wrk'
        if cr.datestr == None:
            self.getDateString(cr)
        datestr = cr.datestr
        # day = str(cr.date.ymdhms['day'])
        # day2 = str(cr.date.ymdhms['day'] + 1)
        # month = str(cr.date.ymdhms['month'])
        # if cr.date.ymdhms['day'] < 10:
        #     day = '0' + str(cr.date.ymdhms['day'])
        #     if cr.date.ymdhms['day'] < 9:
        #         day2 = '0' + str(cr.date.ymdhms['day'] + 1)
        # if cr.date.ymdhms['month'] < 10:
        #     month = '0' + str(cr.date.ymdhms['month'])

        # datestr = str(cr.date.ymdhms['year']) + '-' + month + '-' + day + '+' + day2
        print(datestr)
        self.mkwrk(cr)
        # mk wrk folder, allow for it to already exist
        # os.makedirs(wrkdir / datestr, exist_ok = True)
        # os.makedirs(wrkdir / datestr / 'aligned', exist_ok = True)
        # os.makedirs(wrkdir / datestr / 'calibrated', exist_ok = True)
        # os.makedirs(wrkdir / datestr / 'uncalibrated', exist_ok = True)
        if filters == None:
            filters = cr.filters.keys()
        
        for filter in filters:
            fildat = cr.data[cr.filters[filter]]
            if (fildat.target == None):
                fplate = str(int(cr.date.mjd)) + '-' + filter + '_'
            else:
                fplate = str(fildat.target.name) + '-' + filter + '_' 
            if (fildat.calibrated == True) and (fildat.aligned == True): 
                wrdir = wrkdir / datestr / 'aligned'
                if len(filters) != 1:
                    os.makedirs(wrkdir / datestr / 'aligned' / filter, exist_ok = True)
                    wrdir = wrdir / filter
                fsub = '_ca'
            elif (fildat.calibrated == True):
                wrdir = wrkdir / datestr / 'calibrated'
                if len(filters) != 1:
                    os.makedirs(wrkdir / datestr / 'calibrated' / filter, exist_ok = True)
                    wrdir = wrdir / filter
                fsub = '_c'
            else: 
                wrdir = wrkdir / datestr / 'uncalibrated'
                if len(filters) != 1:
                    os.makedirs(wrkdir / datestr / 'uncalibrated' / filter, exist_ok = True)
                    wrdir = wrdir / filter
                fsub = ''

                
            for p in range(len(fildat.data)):
                image = fildat.data[p]
                fname = fplate + str(p) + fsub + '.fits'
                image.write(wrdir / fname, overwrite = True)

    def saveWCS(self, cr, filters = None):
        wrkdir = self.dordir / 'data' / 'wrk'
        if cr.datestr == None:
            self.getDateString(cr)
        datestr = cr.datestr
        self.mkwrk(cr)
        if filters == None:
            filters = cr.filters.keys()
        
        for filter in filters:
            stack = cr.data[cr.filters[filter]]
            if stack.wcs == None:
                cr.getWCS(filter, self)
            fname = str(filter) + '-solved.fits'
            solved = cr.data[cr.filters[filter]].solved
            solved.write(wrkdir / datestr / 'WCS' / fname, overwrite = True)
            


    # merge header





# we need to pass the calibration files(if none grab the most recent out of 
# dorado by closest mjd)
# if there is no calibration files in dorado send a warning and set 
# calibration files to none

# make a function to force uint16





