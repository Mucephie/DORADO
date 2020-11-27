import dorado
import astroalign as aa


## Organize and prepare your data
# Get a fits file formated date string for the current observing night.
night = dorado.get_night() 
print(night)
# Your data directory
# Set night as directory 
directory = 'C:/Data/' + night

# You can name your target object
target = 'dorado'

# Create a working directory
dorado.mkwrkdir(directory)


# Catalogue input data from the data directory
bias_list, flats_list, lights_list = dorado.checkdir(directory)





## Calibrate your data
# Read the data in
data_series = dorado.get_series(directory, lights_list)
flats_series = dorado.get_series(directory, flats_list)
bias_series = dorado.get_series(directory, bias_list)

# Produce master reduction images
bias = dorado.mastBias(directory, bias_list)
flat = dorado.mastFlat(directory, flats_list, bias)

# Calibrate data series
series = dorado.reduce_series(directory, lights_list, flat, bias, target)

# Write calibrated data
dorado.write_series(directory, series, target)





## Plate solve your data series
# Pick image to solve and store its data
data = series[0]
# Set up a file string for dorado to locate image
expath = '/wrk/calibrated/' 
im = target +  '-0_c.fit'
image_file_path = dorado.file_string(directory + expath, im)
# Pass the data and file string to dorado
solved = dorado.plate_solve(data, image_file_path, write_fits = True, write_name = target + '_solved')




# Align each image with the plate solved image
aa_series = []
for image in series:
    img_r, footprint_r = aa.register(image.data, solved.data)
    aa_series.append(img_r)

dorado.write_series(directory + '/wrk/aligned/', aa_series, target + '_a')
# You now have a fully calibrated and aligned data set with WCS data









