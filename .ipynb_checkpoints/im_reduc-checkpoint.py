import draco
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table



file_path = 'C:/Users/mucep/Offline/Draco-test/'
im_prefix = 'Draco-'
c46P_prefix = '46P_'
BIAS_suffix = '.BIAS'
FLAT_prefix = 'FLAT_'
FLAT_suffix = '.FLAT'
filters = ['BLUE', 'CLEAR', 'EIGHT', 'GREEN', 'LUMINANCE', 'RED', 'SEVEN', 'SIX']

# master bias made by DRACO
master_bias = draco.mast_reduc(file_path + 'Bias/', im_prefix, BIAS_suffix, 10, False)
#fits.writeto('master_bias.fits', master_bias)

reduc_filters =  ['BLUE', 'CLEAR', 'LUMINANCE', 'SEVEN', 'SIX']
size = master_bias.shape
im_count = 5
master_i_Flat = np.zeros((size[0], size[1], len(reduc_filters)))

for i in range(len(reduc_filters)):
    master_i_Flat[:, :, i] = draco.mast_reduc(file_path + 'Flats/' + reduc_filters[i] + '/', im_prefix + reduc_filters[i] + '-', FLAT_suffix, 10, False)
    master_i_Flat[:, :, i] = draco.imarith(master_i_Flat[:, :, i], '-', master_bias)
    master_i_Flat[:, :, i] = draco.norm_flat(master_i_Flat[:, :, i])
    #fits.writeto('master_' + reduc_filters[i] + '_Flat.fits', master_i_Flat[i])

series_n = np.zeros((size[0], size[1], im_count, len(reduc_filters)))
for n in  range(len(reduc_filters)):
    series_n[:, :, :, n] = draco.get_series(file_path + 'Lights/', c46P_prefix + reduc_filters[n] + '-', '', 5, False)
    series_n[:, :, :, n] = draco.series_arith(series_n[:, :, :, n],'-', master_bias, 5)
    series_n[:, :, :, n] = draco.series_arith(series_n[:, :, :, n],'/', master_i_Flat[:, :, n], 5)
    sky = draco.sky_est(series_n[:, :, :, n])
    series_n[:, :, :, n] = draco.series_arith(series_n[:, :, :, n],'-', sky, 5)
    #draco.series_write(series_n[:, :, :, n], 'c46P-' + reduc_filters[n] + '-',5)

draco.plt_fits(series_n[:, :, 3, 3])

## Data reduction checklist:

# produce master bias

# produce master flats

# remove bias from images/lights

# remove bias from flats

# normalize flats

# divide bias corrected lights by normalized flats

# align fully corrected images/lights

# combine aligned images/lights

# estimate and remove sky values

# maybe make a median filter for comet46P and unsharp masked image to expose nuclei