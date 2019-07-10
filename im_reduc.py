
# # IRAF made master bias image
# master_bias_header, master_bias_data = get_im(file_path, 'MasterBias')
# plt_fits(master_bias_data)

# master bias made by DRACO
master_bias = mast_reduc(file_path + 'Bias/', im_prefix, BIAS_suffix, 10)
#fits.writeto('master_bias.fits', master_bias)
#plt_fits(master_bias)

# [BLUE, CLEAR, LUMINANCE, SEVEN, SIX]

# master BLUE Flat made by DRACO
master_BLUE_Flat = mast_reduc(file_path + 'Flats/BLUE/', im_prefix + 'BLUE-', FLAT_suffix, 10)
master_BLUE_Flat = imarith(master_BLUE_Flat, '-', master_bias)
master_BLUE_Flat = norm_flat(master_BLUE_Flat)
#fits.writeto('master_BLUE_Flat.fits', master_BLUE_Flat)
#plt_fits(master_BLUE_Flat)
# master CLEAR Flat made by DRACO
master_CLEAR_Flat = mast_reduc(file_path + 'Flats/CLEAR/', im_prefix + 'CLEAR-', FLAT_suffix, 10)
master_CLEAR_Flat = imarith(master_CLEAR_Flat, '-', master_bias)
master_CLEAR_Flat = norm_flat(master_CLEAR_Flat)
#fits.writeto('master_CLEAR_Flat.fits', master_CLEAR_Flat)
#plt_fits(master_CLEAR_Flat)
# master LUMINANCE Flat made by DRACO
master_LUM_Flat = mast_reduc(file_path + 'Flats/LUMINANCE/', im_prefix + 'LUMINANCE-', FLAT_suffix, 10)
master_LUM_Flat = imarith(master_LUM_Flat, '-', master_bias)
master_LUM_Flat = norm_flat(master_LUM_Flat)
#fits.writeto('master_LUM_Flat.fits', master_LUM_Flat)
#plt_fits(master_LUM_Flat)
# master SEVEN Flat made by DRACO
master_SEVEN_Flat = mast_reduc(file_path + 'Flats/SEVEN/', im_prefix + 'SEVEN-', FLAT_suffix, 10)
master_SEVEN_Flat = imarith(master_SEVEN_Flat, '-', master_bias)
master_SEVEN_Flat = norm_flat(master_SEVEN_Flat)
#fits.writeto('master_SEVEN_Flat.fits', master_SEVEN_Flat)
#plt_fits(master_SEVEN_Flat)
# master SIX Flat made by DRACO
master_SIX_Flat = mast_reduc(file_path + 'Flats/SIX/', im_prefix + 'SIX-', FLAT_suffix, 10)
master_SIX_Flat = imarith(master_SIX_Flat, '-', master_bias)
master_SIX_Flat = norm_flat(master_SIX_Flat)
#fits.writeto('master_SIX_Flat.fits', master_SIX_Flat)
#plt_fits(master_SIX_Flat)




# master BLUE light made by DRACO
series_BLUE = get_series(file_path + 'Lights/', c46P_prefix + 'BLUE-', '', 5)
series_BLUE = series_arith(series_BLUE,'-', master_bias, 5)
series_BLUE = series_arith(series_BLUE,'/', master_BLUE_Flat, 5)
series_BLUE = series_arith(series_BLUE,'-', 54, 5)
#series_write(series_BLUE, 'c46P-BLUE-SKY-',5)
plt_fits(series_BLUE[:, :, 3])
# master CLEAR light made by DRACO
series_CLEAR = get_series(file_path + 'Lights/', c46P_prefix + 'CLEAR-', '', 5)
series_CLEAR = series_arith(series_CLEAR,'-', master_bias, 5)
series_CLEAR = series_arith(series_CLEAR,'/', master_CLEAR_Flat, 5)
#series_write(series_CLEAR, 'c46P-CLEAR-',5)
# master LUMINANCE light made by DRACO
series_LUM = get_series(file_path + 'Lights/', c46P_prefix + 'LUM-', '', 5)
series_LUM = series_arith(series_LUM,'-', master_bias, 5)
series_LUM = series_arith(series_LUM,'/', master_LUM_Flat, 5)
#series_write(series_LUM, 'c46P-LUM-',5)
# master SEVEN light made by DRACO
series_SEVEN = get_series(file_path + 'Lights/', c46P_prefix + 'SEVEN-', '', 5)
series_SEVEN = series_arith(series_SEVEN,'-', master_bias, 5)
series_SEVEN = series_arith(series_SEVEN,'/', master_SEVEN_Flat, 5)
#series_write(series_SEVEN, 'c46P-SEVEN-',5)
# master SIX light made by DRACO
series_SIX = get_series(file_path + 'Lights/', c46P_prefix + 'SIX-', '', 5)
series_SIX = series_arith(series_SIX,'-', master_bias, 5)
series_SIX = series_arith(series_SIX,'/', master_SIX_Flat, 5)
#series_write(series_SIX, 'c46P-SIX-',5)

#plt_fits(master_BLUE)




# remove bias from images

# remove bias from flats

# normalize flats

# divide bias corrected lights by normalized flats

# align fully corrected images

# combine aligned images

# estimate and remove sky values

# maybe make a median filter for comet46P and unsharp masked image to expose nuclei