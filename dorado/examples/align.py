import dorado
import astroalign as aa
import datetime

## One might be interested in how long this all takes
START_DATE_TIME = datetime.datetime.now()
print('\nStarting time: ', START_DATE_TIME)
print('\nFinding data...')



## Organize and prepare your data
# Get a fits file formated date string for the current observing night.
night = dorado.get_night() 
print('Night: ', night)
# Your data directory
# Set night as directory 
directory = 'C:/Data/' + night

# You can name your target object
target = 'WASP33b'

bias_list, flats_list, lights_list = dorado.checkdir(directory + '/wrk/calibrated/')

# Read the data in
series = dorado.get_series(directory + '/wrk/calibrated/', lights_list[0:200], unit = None)

CAL_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', CAL_END_DATE_TIME)
print("Time elapsed for sheet: ", (CAL_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')

## Plate solve your data series
# Pick image to solve and store its data
psimg = series[0]
# Set up a file string for dorado to locate image
expath = '/wrk/calibrated/' 
im = target +  '-0_c.fit'
image_file_path = dorado.file_string(directory + expath, im)
# Pass the data and file string to dorado
solved = dorado.plate_solve(psimg, image_file_path, write_fits = True, write_name = target + '_solved')
PLATE_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', PLATE_END_DATE_TIME)
print("Time elapsed for sheet: ", (PLATE_END_DATE_TIME - CAL_END_DATE_TIME))
print("Time elapsed so far: ", (PLATE_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')





# Align each image with the plate solved image
aa_series = []
for image in series:
    img_r, footprint_r = aa.register(image.data, solved.data)
    aa_series.append(img_r)

dorado.write_series(directory + '/wrk/aligned/', aa_series, target, '_a')
# You now have a fully calibrated and aligned data set with WCS data
AL_END_DATE_TIME = datetime.datetime.now()
print('\nFinished time: ', AL_END_DATE_TIME)
print("Time elapsed for sheet: ", (AL_END_DATE_TIME - PLATE_END_DATE_TIME))
print("Time elapsed so far: ", (AL_END_DATE_TIME - START_DATE_TIME))
print('\n\n\n\n')
