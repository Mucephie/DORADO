#ts = QTable([times, exptimes, x, y, ray, decx, flux, fluxunc], names=('time', 'exptime', 'x', 'y', 'ra', 'dec', 'flux', 'flux_unc'), meta={'name': filter})

import numpy as np
from astropy.time import Time
from astropy.table import QTable, Table
from astropy.io import fits

__all__ = ['timeSeries']



class timeSeries:
    '''
        The timeSeries

        Attributes
        ----------

        name: str
            name of target in string format.

    '''
    def __init__(self, times, flux, exptimes = [], x = [], y = [], ra = [], dec = [], flux_unc = [], apsum = [], apsum_unc = [], fit_times = [], fit_flux = [], toml = [], OmC = [], cycle = []):  
        self.times = times
        self.flux = flux

        self.exptimes = exptimes

        self.x = x
        self.y = y
        self.ra = ra
        self.dec = dec

        self.flux_unc = flux_unc
        self.apsum = apsum
        self.apsum_unc = apsum_unc

        self.fit_times = fit_times
        self.fit_flux = fit_flux

        self.toml = toml
        self.OmC = OmC
        self.cycle = cycle

        # self.symbo = None # make a symbolic expression to represent the curve analytically.

    def toTable(self, name):

        self.table = QTable([self.times, self.flux], names=('time','flux'), meta={'name': name})

        colnom = ['flux_unc', 'exptime', 'apsum', 'apsum_unc', 'x', 'y', 'ra', 'dec']
        cols = [self.flux_unc, self.exptimes, self.apsum, self.apsum_unc, self.x, self.y, self.ra, self.dec]

        for col in range(len(colnom)):
            if cols[col] != []:
                try:
                    self.table[colnom[col]] = cols[col]
                except:
                    print('Error merging', colnom[col], ' with table.')