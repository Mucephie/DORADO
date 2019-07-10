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



# # IRAF made master bias image
# master_bias_header, master_bias_data = get_im(file_path, 'MasterBias')
# plt_fits(master_bias_data)

# master bias made by DRACO
master_bias = draco.mast_reduc(file_path + 'Bias/', im_prefix, BIAS_suffix, 10)
#fits.writeto('master_bias.fits', master_bias)
#plt_fits(master_bias)

reduc_filters =  ['BLUE', 'CLEAR', 'LUMINANCE', 'SEVEN', 'SIX']

# master BLUE Flat made by DRACO
master_BLUE_Flat = draco.mast_reduc(file_path + 'Flats/BLUE/', im_prefix + 'BLUE-', FLAT_suffix, 10)
master_BLUE_Flat = draco.imarith(master_BLUE_Flat, '-', master_bias)
master_BLUE_Flat = draco.norm_flat(master_BLUE_Flat)
#fits.writeto('master_BLUE_Flat.fits', master_BLUE_Flat)

# master CLEAR Flat made by DRACO
master_CLEAR_Flat = draco.mast_reduc(file_path + 'Flats/CLEAR/', im_prefix + 'CLEAR-', FLAT_suffix, 10)
master_CLEAR_Flat = draco.imarith(master_CLEAR_Flat, '-', master_bias)
master_CLEAR_Flat = draco.norm_flat(master_CLEAR_Flat)
#fits.writeto('master_CLEAR_Flat.fits', master_CLEAR_Flat)

# master LUMINANCE Flat made by DRACO
master_LUM_Flat = draco.mast_reduc(file_path + 'Flats/LUMINANCE/', im_prefix + 'LUMINANCE-', FLAT_suffix, 10)
master_LUM_Flat = draco.imarith(master_LUM_Flat, '-', master_bias)
master_LUM_Flat = draco.norm_flat(master_LUM_Flat)
#fits.writeto('master_LUM_Flat.fits', master_LUM_Flat)

# master SEVEN Flat made by DRACO
master_SEVEN_Flat = draco.mast_reduc(file_path + 'Flats/SEVEN/', im_prefix + 'SEVEN-', FLAT_suffix, 10)
master_SEVEN_Flat = draco.imarith(master_SEVEN_Flat, '-', master_bias)
master_SEVEN_Flat = draco.norm_flat(master_SEVEN_Flat)
#fits.writeto('master_SEVEN_Flat.fits', master_SEVEN_Flat)

# master SIX Flat made by DRACO
master_SIX_Flat = draco.mast_reduc(file_path + 'Flats/SIX/', im_prefix + 'SIX-', FLAT_suffix, 10)
master_SIX_Flat = draco.imarith(master_SIX_Flat, '-', master_bias)
master_SIX_Flat = draco.norm_flat(master_SIX_Flat)
#fits.writeto('master_SIX_Flat.fits', master_SIX_Flat)





# master BLUE light made by DRACO
series_BLUE = draco.get_series(file_path + 'Lights/', c46P_prefix + 'BLUE-', '', 5)
series_BLUE = draco.series_arith(series_BLUE,'-', master_bias, 5)
series_BLUE = draco.series_arith(series_BLUE,'/', master_BLUE_Flat, 5)
series_BLUE = draco.series_arith(series_BLUE,'-', 54, 5)
#draco.series_write(series_BLUE, 'c46P-BLUE-SKY-',5)
draco.plt_fits(series_BLUE[:, :, 3])
# master CLEAR light made by DRACO
series_CLEAR = draco.get_series(file_path + 'Lights/', c46P_prefix + 'CLEAR-', '', 5)
series_CLEAR = draco.series_arith(series_CLEAR,'-', master_bias, 5)
series_CLEAR = draco.series_arith(series_CLEAR,'/', master_CLEAR_Flat, 5)
#draco.series_write(series_CLEAR, 'c46P-CLEAR-',5)
# master LUMINANCE light made by DRACO
series_LUM = draco.get_series(file_path + 'Lights/', c46P_prefix + 'LUM-', '', 5)
series_LUM = draco.series_arith(series_LUM,'-', master_bias, 5)
series_LUM = draco.series_arith(series_LUM,'/', master_LUM_Flat, 5)
#draco.series_write(series_LUM, 'c46P-LUM-',5)
# master SEVEN light made by DRACO
series_SEVEN = draco.get_series(file_path + 'Lights/', c46P_prefix + 'SEVEN-', '', 5)
series_SEVEN = draco.series_arith(series_SEVEN,'-', master_bias, 5)
series_SEVEN = draco.series_arith(series_SEVEN,'/', master_SEVEN_Flat, 5)
#draco.draco.series_write(series_SEVEN, 'c46P-SEVEN-',5)
# master SIX light made by DRACO
series_SIX = draco.get_series(file_path + 'Lights/', c46P_prefix + 'SIX-', '', 5)
series_SIX = draco.series_arith(series_SIX,'-', master_bias, 5)
series_SIX = draco.series_arith(series_SIX,'/', master_SIX_Flat, 5)
#draco.series_write(series_SIX, 'c46P-SIX-',5)

#draco.plt_fits(master_BLUE)




# remove bias from images

# remove bias from flats

# normalize flats

# divide bias corrected lights by normalized flats

# align fully corrected images

# combine aligned images

# estimate and remove sky values

# maybe make a median filter for comet46P and unsharp masked image to expose nuclei