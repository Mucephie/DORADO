import warnings
warnings.filterwarnings('ignore')
import os
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord as acoord
import numpy as np
import astropy.units as un

class zellars:
    def __init__(self, name):
        # get it because zellars is the canadian target?
        self.name = name
        s = Simbad()
        r = s.query_object(self.name)
        self.filters = {}
        self.ts = []
        # r.pprint()
        # print(r.colnames)
        self.coords = acoord(ra = r['RA'], dec = r['DEC'], unit = (un.hourangle, un.deg), frame = 'icrs')
        # xy coordinates
        
    def calcmag(self, filter):
        flux  = self.ts[self.filters[filter]]['flux']
        flux_unc = self.ts[self.filters[filter]]['flux_unc']
        magnitudes = -2.5 * np.log10(flux / 15)
        mag_unc = flux_unc / ((flux / 15) * 2.30258509)
        self.ts[self.filters[filter]]['mag'] = magnitudes
        self.ts[self.filters[filter]]['mag_unc'] = mag_unc

    def record(self, clippy, cr, saveType = 'fits'):
        wrkdir = clippy.dordir / 'data' / 'wrk'
        if cr.datestr == None:
            cr.datestr = clippy.getDateString(cr)
        datestr = cr.datestr
        wrdir = wrkdir / datestr
        clippy.mkwrk(cr)
        for fi in self.filters.keys():
            wrts = self.ts[self.filters[fi]]
            fname = str(self.name) + '_' + str(fi) + '-' + str(int(cr.date.mjd)) + '.' + saveType
            wrts.write(wrdir / fname, overwrite = True)
        
    def export(self, clippy, cr, objectClass = None):
        tardir = clippy.dordir / 'targets'
        if objectClass:
            os.makedirs(tardir / objectClass, exist_ok = True)
            tardir = tardir / objectClass
        
