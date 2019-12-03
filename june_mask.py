#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 14:32:08 2019

@author: benmendez
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from photutils import aperture_photometry
from photutils import RectangularAperture
from photutils import RectangularAnnulus

#aperture = RectangularAperture((6,5), 3,3)
#
#annulus_aperture = RectangularAnnulus((6,5),4,5,5)

#aperture = RectangularAperture((6,5), 4,4)
#
#annulus_aperture = RectangularAnnulus((6,5),4,5,5)
#
#apers = [aperture , annulus_aperture]

#aperture = CircularAperture(positions, r=3)

#masks = aperture.to_mask(method='center')
#mask = masks[0]
#image = mask.to_image(shape=((11, 11)))
#data_cutout = mask.cutout(calibrated_fluxes[a,:,:])
#data_cutout_aper = masks.multiply(calibrated_fluxes[a,:,:])

for a in range(calibrated_fluxes.shape[0]):
    
    #Setting up the apertures around KELT-16
    aperture = RectangularAperture((6,5), 3,3)
    annulus_aperture = RectangularAnnulus((6,5), 4, 5, 5)
    apers = [aperture , annulus_aperture]
    
    #Applying TESS data through the apertures to determine flux through apertures
    phot_table = aperture_photometry(calibrated_fluxes[a,:,:],apers)
    
    #Iterating over all cadences
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    
    #Calculating background noise through aperture
    bkg_mean = phot_table['aperture_sum_1']/annulus_aperture.area
    bkg_sum = bkg_mean * aperture.area
    
    #Subtracting background from main aperture
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum
    
    #Appending final flux value into array to be used for plotting light curve
    aperature_sums_a.append(final_sum)
    
    
    phot_table['residual_aperture_sum'].info.format = '%.8g'  # for consistent table output
    ap_mean = phot_table['aperture_sum_0']/aperture.area
    
    #print(phot_table['residual_aperture_sum'])  
    
   # phot_table['residual_aperture_sum'] = final_sum
    
  #  aperature_sums_b.append(2.5 * np.log10(converted_aperture_sum))
  #  aperature_sums_b.append(converted_aperture_sum)
 #   aperature_sums_a.append(final_sum)
    #mean_1.append(bkg_sum)
