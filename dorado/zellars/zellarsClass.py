import warnings
warnings.filterwarnings('ignore')
import os
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord as acoord
import numpy as np
import astropy.units as un

__all__ = ['Zellars']

class Zellars:
    '''
        The Zellars class represents an astronomical target and handles the targets relevent attributes.
        Unpassed target parameters will be gathered via astroquery (SIMBAD).

        The name is a reference to the now defunct Canadian store chain, Zellars, which many Canadians
        saw as the Canadian version of the American store chain, Target.

        Attributes
        ----------

        name: str
            name of target in string format.
    '''
    def __init__(self, name):
        self.name = name
        # Call Simbad for relevent data
        s = Simbad()
        r = s.query_object(self.name)
        self.coords = acoord(ra = r['RA'], dec = r['DEC'], unit = (un.hourangle, un.deg), frame = 'icrs')
        
        self.filters = {}
        self.ts = []

        
    def calcmag(self, filter):
        '''
        calcmag converts the targets flux and associated uncertainty into an 
        instrumental magnitude and uncertainty.

        Parameters
        ----------
        filter: str
            String representation of the relevent filter.

        Returns
        -------
        self.ts[filter]['mag'], self.ts[filter]['mag_unc']

        '''
        flux  = self.ts[self.filters[filter]]['flux']
        flux_unc = self.ts[self.filters[filter]]['flux_unc']
        magnitudes = -2.5 * np.log10(flux / 15)
        mag_unc = flux_unc / ((flux / 15) * 2.30258509)
        self.ts[self.filters[filter]]['mag'] = magnitudes
        self.ts[self.filters[filter]]['mag_unc'] = mag_unc

    def record(self, clippy, cr, saveType = 'fits'):
        '''
        record writes each filters timeseries to the dorado working data directory
        for the relevent date and target.

        Parameters
        ----------
        clippy: Clippy instance
            The active instance of clippy to handle writing the file

        cr: Ceres instance
            The relevent instance of Ceres for which the timeseries was derived. The save location 
            will be the working directory for this instance.
            
        saveType: str
            String representation of the file extension to save to. Default is 'fits'. See 
            'astropy.table - write()' for acceptable values.

        Returns
        -------
        None
        '''
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
        
    def export(self, clippy, objectClass = None):
        '''
        export will record the Zellars object into the dorado targets directory for future use 
        and reference. This function is not implemented yet.

        Parameters
        ----------
        clippy: Clippy instance
            The active instance of clippy to handle writing the target file
        objectClass: str
            Class of object to save the target under. 

        Notes
        -----
        Examples of object classes may be: 'star', 'galaxy', 'exoplanet', 'minor planet', 'satellite', 
        'white dwarf', 'nebula', 'messier_object', 'O_star', 'binary', 'globular_cluster', 'open_cluster', 
        'galaxy_cluster', 'quasar', 'AGN'. 
        
        Users can craft their own object naming schemes.
        '''
        tardir = clippy.dordir / 'data/targets'
        if objectClass:
            os.makedirs(tardir / objectClass, exist_ok = True)
            tardir = tardir / objectClass
        
