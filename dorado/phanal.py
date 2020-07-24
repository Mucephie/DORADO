##  imports
import numpy as np
from scipy.ndimage import median_filter
from skimage.feature import blob_dog, blob_log
import skimage.exposure as skie

__all__ = ['starSeeker']

def starSeeker(data):
        mf = median_filter(data, size= 15)
        datamf = data - mf
        limg = np.arcsinh(datamf) #datamf
        limg = limg / limg.max()
        low = np.percentile(limg, 0.2)
        high = np.percentile(limg, 99.1)
        opt_img  = skie.exposure.rescale_intensity(limg, in_range=(low,high))

        stars =  blob_log(opt_img, max_sigma=100, min_sigma = 5, num_sigma=10, threshold=.2)
        # stars =  blob_dog(opt_img, max_sigma=30, threshold=.2)

        # Compute radii in the 3rd column.
        stars[:, 2] = stars[:, 2] * np.sqrt(2)

        y2, x2, r = stars[:, 0], stars[:, 1], stars[:, 2]


        print('\nStars found: ', len(r))

        return x2, y2, r, opt_img

# def apphot
# def diff_phot
# def lc_clean
# def waldo