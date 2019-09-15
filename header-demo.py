from astropy.io import fits
import numpy as np
from astropy.table import Table

from astropy.nddata import CCDData
from ccdproc import Combiner
import ccdproc
import astropy.units as u
from astropy.utils.misc import isiterable

import draco

directory = 'C:/Users/mucep/Offline/2019-08-14+15'
star = 'YZBoo'
bias_list, flats_list, lights_list = draco.checkdir(directory)
bias, bias_headers = draco.get_series(directory, bias_list)
# master_bias = draco.mast_reduc(bias)

flats, flat_headers = draco.get_series(directory, flats_list)
print((flat_headers[0][29])) #._keyword_indices)