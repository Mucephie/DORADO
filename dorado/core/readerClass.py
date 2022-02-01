import warnings
warnings.filterwarnings('ignore')

from ..ceres import Ceres
from ..stack import Stack
from ..core.coreClass import *

import ccdprocx

from astropy.nddata.ccddata import CCDData
from astropy.time import Time
from astropy.io import fits

import lightkurve as lk
from lightkurve.lightcurve import TessLightCurve as tlc

import os

__all__ = ['aico_reader'] #, 'tess_reader'

# TODO Make it so that each class isnt loaded unless chosen (core class takes string and then imports correct class)
class aico_reader:
    def __init__(self):
        self.temp = None

    def dirscan(self, dirarray):
        # TODO handle multi filter searching 
        path = Dorado.dordir
        for dir in dirarray:
            path = path / dir
        files, directories = Dorado.diread(path)

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
                    hdu = CCDData.read(biasl[i].path, unit = Dorado.unit)
                    bias.append(hdu)
                # print('Bias searched.')
                flatsl = []
                for strflat in flatsstr:
                    for s in files:
                        if strflat in s.name:
                            flatsl.append(s)
                flats = []
                for i in flatsl:
                    hdu = CCDData.read(i.path, unit = Dorado.unit)
                    flats.append(hdu)
                # print('flats searched.')
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
                    hdu = CCDData.read(i.path, unit = Dorado.unit)
                    ## Trying to enforce 16bit data instead of 32 or 64
                    hdu.data = hdu.data.astype('uint16')
                    lights.append(hdu)
                # print('lights searched.')
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
            biasl, _ = Dorado.diread(biasdir[0])
            bias = []
            for i in biasl:
                hdu = CCDData.read(i.path, unit = Dorado.unit)
                bias.append(hdu)


            for ldir in lightsdir:
                files, directories = Dorado.diread(ldir)
                if len(directories) == 0:
                    if len(files) == 0:
                        raise Exception('No viable light data found')
                    else:
                        print('Single directory lights organization format detected.')
                        # print('Reading files.')
                        lights = []
                        for i in files:
                            hdu = CCDData.read(i.path, unit = Dorado.unit)
                            ## Trying to enforce 16bit data instead of 32 or 64
                            hdu.data = hdu.data.astype('uint16')
                            lights.append(hdu)
                        # filter = ImageFileCollection(ldir).values('filter', unique = True)
                        # lights = [filter, lightsarr]

                        
                elif len(files) == 0:
                    print('Multi directory lights organization format detected.')
                    # print('Reading directories.')
                    lights = []



            for fdir in flatsdir:
                files, directories = Dorado.diread(fdir)
                if len(directories) == 0:
                    if len(files) == 0:
                        raise Exception('No viable light data found')
                    else:
                        print('Single directory flats organization format detected.')
                        # print('Reading files.')
                        flats = []
                        for i in files:
                            hdu = CCDData.read(i.path, unit = Dorado.unit)
                            flats.append(hdu)

                elif len(files) == 0:
                    print('Multi directory flats organization format detected.')
                    # print('Reading directories.')
                    flats = []

           
            return bias, flats, lights
        
    def mkceres(self,  date, sub = 'raw', target = None, calibrated = False, aligned = False):
            # TODO needs ability to handle multifilter directories
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
            ## TODO :: get missing calibraation files
            ## if data aligned or calibrated, cal frames arent needed

            ## TODO :: look into UTC wrecking stuff
            if len(biasIFC) == 0:
                cere = Ceres(time = Time(lights[0].header['DATE-OBS'], format='fits'))
                ## TODO get date string is in filer
                # self.getDateString(cere)
            else:
                bias = self.mkBias(biasIFC)
                cere = Ceres(bias = bias, time = Time(lights[0].header['DATE-OBS'], format='fits'))
                # self.getDateString(cere)

            # for f in filter_data:
            # cere.flats[flat.header['filter']] = flat

            ## TODO :: multifilter fun
            if len(flats) == 0:
                 cere.add_stack(Stack(lights, calibrated = calibrated, aligned = aligned, target = target))
            else:
                flat = self.mkFlat(flats)
                cere.add_stack(Stack(lights, flat = flat, calibrated = calibrated, aligned = aligned, target = target))
            

            return cere
        
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
            # Allow RGB data
            c = ccdprocx.Combiner(flats)
            c.sigma_clipping()
            flat = c.median_combine()
            # , method = 'average',
            #                     sigma_clip = True, sigma_clip_low_thresh = 5, sigma_clip_high_thresh = 5,
            #                     sigma_clip_func = np.ma.median, sigma_clip_dev_func = mad_std, unit = self.unit)
            flat.header['stacked'] = True
            flat.header['numsubs'] = len(flats)
            flat.header['DATE-OBS'] = flats[0].header['DATE-OBS']
            flat.header['filter'] = flats[0].header['filter']
            ## TODO :: There is probably more missing keywords in the combined header, where's Waldo...

            date = Time(flat.header['DATE-OBS'], format='fits').mjd
            flat.data = flat.data.astype('uint16') 

            filt = flat.header['filter']
            ## TODO :: standardize filter names and include telescope profiles
            fname = str(int(date)) + '_' + str(filt) + '_flat.fits'
            flatdir = Dorado.dordir / 'data' / 'flats' 
            contents = os.scandir(path = flatdir)
            save = True
            for entry in contents:
                if fname in entry.name:
                    save = False
            if save:
                print('Saving Flat for later use')
                flat.write(flatdir / fname)
            else:
                print('Flat for filter and date already saved.')

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
            # Allow RGB data

            bias = ccdprocx.combine(biasIFC, method = 'average', unit = Dorado.unit)
            bias.meta['stacked'] = True
            bias.header['numsubs'] = len(biasIFC)
            date = Time(bias.header['DATE-OBS'], format='fits').mjd
            bias.data = bias.data.astype('uint16') 
            fname = str(int(date)) + '_Bias.fits'
            biasdir = Dorado.dordir / 'data' / 'bias' 
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
        


# class tess_reader:
#     def __init__(self):
#         Dorado.init_dir(tess = True)
#         # check if a name query is needed to retrieve data
#         # initialize 

#     def mkceres(self, tid, sub = 'raw', target = None):
#         fits_file = self.dirscan()
#         tpf = lk.open(fits_file)

#         with fits.open(fits_file, mode="readonly") as hdulist:
#             tess_bjds = hdulist[1].data['TIME'] # Should I use 'TIMECORR'?
#             tess_bjds_corr = hdulist[1].data['TIMECORR'] 
#             raw_counts = hdulist[1].data['RAW_CNTS']
#             calibrated_fluxes = hdulist[1].data['FLUX'] # Should I use raw or calibrated?
#             flux_err = hdulist[1].data['FLUX_ERR']
#             tid = hdulist[0].header['TICID']
#             TSTART = hdulist[0].header['TSTART']
#             TSTOP = hdulist[0].header['TSTOP']
#             tidstr = 'TIC ' + str(tid)
        
#         cere = Ceres()
#         # need to make a stack, and reconfigure the ceres class to be more generic in data storage. maybe make the dorphot function
#         # a wrapper, or make the ability to produce a target pixel file of aico data.
#         return cere


#     def dirscan(self, Dorado, dirarray):
#         path = Dorado.dordir
#         for dir in dirarray:
#             path = path / dir

#         files, directories = Dorado.diread(path)

#         if len(files) == 0:
#             print('No files found...')
#         elif len(files) > 1:
#             for file in files:
#                 if '_tp.fits' in file:
#                     tp_file = file    
#         elif len(files) == 1:
#             tp_file = files[0]

#         return tp_file

