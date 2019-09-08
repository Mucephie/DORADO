from astropy.io import fits
import numpy as np
# import matplotlib.pyplot as plt
from astropy.table import Table
# from matplotlib import cm
# from astropy.visualization import astropy_mpl_style
# from astropy.utils.data import get_pkg_data_filename
# from matplotlib.colors import LogNorm
# import matplotlib.colors as colors

from astroscrappy import detect_cosmics
from astropy.nddata import CCDData
from ccdproc import Combiner
import ccdproc
import astropy.units as u
from astropy.utils.misc import isiterable

# from scipy.stats import mode