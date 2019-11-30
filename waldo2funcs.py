import warnings
warnings.filterwarnings('ignore')

## file imports
import os
from astropy.io import fits 
from astropy.nddata import CCDData
from astropy.table import Table
from astropy.nddata import NDData
# from pathlib import Path

## Processing imports
import numpy as np
import ccdproc
from astropy.stats import mad_std
import astroalign
from photutils.psf import extract_stars

## photometry imports
from astropy import units as u

## plotting imports
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm

## other imports
import datetime
import dracoOP2

def dypeg_waldo(data):
        x, y, r, opt_img = dracoOP2.starSeeker2(data)
        if (len(x) < 3):
                return 0, 0
        size = 100
        hsize = (size - 1) / 2

        mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y < (data.shape[0] -1 - hsize)))

        stars_tbl = Table()
        stars_tbl['x'] = x[mask]
        stars_tbl['y'] = y[mask]
        nddata = NDData(data=data)
        starsq = extract_stars(nddata, stars_tbl, size=50)

        l = []
        for u in range(len(x[mask])):
                l.append(np.sum(starsq[u]))

        stars = np.zeros((len(x[mask]), 4))
        stars[:, 0] = x[mask]
        stars[:, 1] = y[mask]
        stars[:, 2] = r[mask]
        stars[:, 3] = l


        stars = stars[stars[:,3].argsort()]
        stars = np.flip(stars, axis = 0)
        if (len(stars) < 2):
                return 0, 0
        star5 = stars[(0,1),:]

        starT = Table()
        starT['x'] = star5[:,0]
        starT['y'] = star5[:,1]
        starT['r'] = star5[:,2]
        starT['l'] = star5[:,3]
        
        return starT, star5

def makePSF(data, stars):
        print('TO-DO')

