# TIC: 236445129

import numpy as np
import matplotlib.pyplot as plt
from lightkurve.lightcurve import LightCurve as LC
#from bls import BLS
from matplotlib import gridspec
import eleanor
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle

from mpl_toolkits.axes_grid1 import make_axes_locatable


star = eleanor.Source(gaia = 1864885215233116032) # TIC: 236445129 234503282 tic = 236445129, sector = 15, 

#
# Full target list for orbit 37, sector 15
# Created Wed Aug 14 21:54:32 EDT 2019 
# Concatenation of GI astero exo bright DDT ppas lists
#
# TICID Camera CCD Tmag RA Dec
# 236444701 1 3 10.98 314.2858 31.5229

print('TIC ID: ', star.tic)
print('Coordinates: ', star.coords)
print('Camera, Chip Location: ', star.camera, ',', star.chip)
print('TESS mag: ', star.tess_mag)


data = eleanor.TargetData(star, height=15, width=15, do_psf=True, do_pca=True)

plt.figure(figsize=(14,6))
plt.imshow(data.post_obj.flux[:,:,0], origin='lower', vmax=200)
plt.title('Postcard', fontsize=18)
plt.show()

star.postcard

fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(14,6))
ax1.imshow(data.tpf[0], origin='lower', vmax=200)
ax1.set_title('TPF', fontsize=20)
ax2.imshow(data.aperture, origin='lower')
ax2.set_title('Best Aperture', fontsize=20)
plt.show()