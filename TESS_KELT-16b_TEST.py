import datetime
import warnings
warnings.filterwarnings('ignore')
#########

from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
from lightkurve.lightcurve import TessLightCurve as tlc
from photutils import aperture_photometry, RectangularAperture, RectangularAnnulus


START_DATE_TIME = datetime.datetime.now()

print('\nStarting time: ', START_DATE_TIME)

# For the purposes of this tutorial, we just know the MAST URL location of the file we want to examine.
fits_file = "tess2019226182529-s0015-0000000236445129-0151-s_tp.fits"
fits.info(fits_file)
print('\n')
fits.getdata(fits_file, ext=1).columns
tpf = lk.open(fits_file)
print('\n')

with fits.open(fits_file, mode="readonly") as hdulist:
    tess_bjds = hdulist[1].data['TIME']
    raw_counts = hdulist[1].data['RAW_CNTS']
    calibrated_fluxes = hdulist[1].data['FLUX']

    flux_err = hdulist[1].data['FLUX_ERR'][0,:,:]

    tid = 'TIC 236445129'

print(type(calibrated_fluxes))
print(calibrated_fluxes.shape)
print('\n')
print(tess_bjds.shape)
print('\n')
print(tid)
print('\n')
print(flux_err.shape)
print('\n')

flx_err = np.mean(flux_err)


with fits.open(fits_file, mode="readonly") as hdulist:
    aperture = hdulist[2].data



# tpf.plot(aperture_mask=tpf.pipeline_mask)

aperture_mask = tpf.create_threshold_mask(threshold=12)

# Plot that aperture
# tpf.plot(aperture_mask=aperture_mask)

# Homebuilt aperture
ap_mask_custom = tpf.pipeline_mask + aperture_mask

# tpf.plot(aperture_mask=ap_mask_custom)



lc = tpf.to_lightcurve(aperture_mask=ap_mask_custom)

lc.errorbar()

flat_lc = lc.flatten(window_length=1001)
#flat_lc.errorbar()
clipped_lc = flat_lc.remove_outliers(sigma=4)
clipped_lc.scatter()
lc_custom = tpf.extract_aperture_photometry(aperture_mask = ap_mask_custom)
lc_custom = lc_custom.remove_outliers(sigma=4) #.flatten(window_length=1001)

lc_custom.scatter()


plt.show()


err = []
aperature_sums = []

for s in range(calibrated_fluxes.shape[0]):
        #Make the apertures
        aperture = RectangularAperture((6,5), 3,3)
        annulus_aperture = RectangularAnnulus((6,5), 4, 5, 5)
        aps = [aperture , annulus_aperture]
        
        
        phot_table = aperture_photometry(calibrated_fluxes[s,:,:], aps)
        # Background subtract
        bkg_sum = (phot_table['aperture_sum_1'] / annulus_aperture.area) * aperture.area
        
        #Subtracting background from main aperture
        phot_table['residual_aperture_sum'] = phot_table['aperture_sum_0'] - bkg_sum
        #phot_table['aperature_magnitude'] = (-2.5 * np.log10(phot_table['residual_aperture_sum']))
        aperature_sums.append(phot_table['residual_aperture_sum'])
 


lc_june = lk.LightCurve(time = tess_bjds, flux = aperature_sums,  flux_err = None, flux_unit = None, time_format=None, time_scale=None, targetid = tid, label = 'KELT-16b', meta=None)
lc_june.scatter()
# lc_june.to_fits( path = "tess2019226182529-s0015-0000000236445129-0151-s_june_lc.fits", overwrite=False, flux_column_name='FLUX')



# # Start figure and axis.
fig, ax = plt.subplots()

# Let's define a title for the figure.
fig.suptitle("KELT 16 b Lightcurve - Sector 15, 2 Min Cadence")
ax.scatter(tess_bjds, aperature_sums)
plt.show()





# ## ----------------------------------------------------------------------------------------
# # # Start figure and axis.
# fig, ax = plt.subplots()

# # Display the calibrated fluxes as an image for the fifth cadence.
# cax = ax.imshow(calibrated_fluxes[40,:,:], cmap=plt.cm.YlGnBu_r, origin="lower")

# # Let's define a title for the figure.
# fig.suptitle("KELT 16 b Calibrated Fluxes - Sector 15, Fifth Cadence")

# # Add a color bar.
# cbar = fig.colorbar(cax)
# plt.show()



# # Start figure and axis.
# fig, ax = plt.subplots()

# # Display the pixels as an image.
# cax = ax.imshow(aperture, cmap=plt.cm.YlGnBu_r, origin="lower")

# # Add a color bar.
# cbar = fig.colorbar(cax)

# # Add a title to the plot.
# fig.suptitle("KELT 16 b Aperture - Sector 15")
# plt.show()



# # Start figure and axis.
# fig, ax = plt.subplots()

# # Display, as an image, the 11x11 table that records the bitmask value of 2 being set.
# cax = ax.imshow(bitmask2_set, cmap=plt.cm.YlGnBu_r, origin="lower")

# # Add a color bar.
# cbar = fig.colorbar(cax)

# # Add a title to the plot.
# fig.suptitle("KELT 16 b Aperture - Sector 15 - Pixels Used In Phot. Ap.")
# plt.show()
END_DATE_TIME = datetime.datetime.now()
print('\nEnding time: ', END_DATE_TIME)
print("Time elapsed: ", (END_DATE_TIME - START_DATE_TIME))