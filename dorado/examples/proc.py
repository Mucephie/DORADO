import dorado
import astroalign
from donuts import Donuts

## Organize and prepare your data
# Get a fits file formated date string for the current observing night.
night = dorado.get_night() 

# Your data directory
directory = 'C:/Data/'

# You can name your target object
target = 'dorado'

# Create a working directory
dorado.mkwrkdir(directory, night)

# Catalogue input data from the data directory
bias_list, flats_list, lights_list = dorado.checkdir(directory, night)

## Calibrate your data
# Read the data in
data_series = dorado.get_series(directory, night, lights_list, unit = None)
flats_series = dorado.get_series(directory, night, flats_list, unit = None)
bias_series = dorado.get_series(directory, night, bias_list, unit = None)

# Produce master reduction images
bias = dorado.mastBias(directory, night, bias_list)
flat = dorado.mastFlat(directory, night, flats_list, bias)

# Calibrate data series
series = dorado.reduce_series(directory, night, lights_list, flat, bias, target)

# Write calibrated data
dorado.write_series(directory, night, series, target)

## Plate solve your data series
# Pick image to solve and store its data
data = series[0].data
# Set up a file string for dorado to locate
expath = '/wrk/calibrated/' 
im = target +  '_0.calibrated.fit'
image_file_path = dorado.file_string(directory, night, im, expath = '')
# Pass the data and file string to dorado
solved = dorado.plate_solve(data, image_file_path, write_fits = True, write_name = '')

# Align each image with the solved image
p, (pos_0, pos_1) = astroalign.find_transform(series[0].data, series[20].data)

for (x1, y1), (x2, y2) in zip(pos_0, pos_1):
    print("({:.2f}, {:.2f}) in source --> ({:.2f}, {:.2f}) in target"
          .format(x1, y1, x2, y2))

# Set up a reference image
refim = lights_list[0]
d = Donuts(refimage = refim, image_ext = 0, subtract_bkg = True, ntiles = 64)
# for each image, compute the x/y translation required to align the images onto the reference image
xp, yp = [], []
for image in lights_list:
    shift_result = d.measure_shift(image)
    x, y = shift_result.x, shift_result.y
    # y = shift_result.y
    xp.append(x.value)
    yp.append(y.value)
    # Also check out shift_result.sky_background
    # print(x, y)