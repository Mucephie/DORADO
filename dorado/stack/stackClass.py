import warnings
warnings.filterwarnings('ignore')
from astropy.time import Time

__all__ = ['Stack']

''' 
'dorado.stack' holds the base Stack class for handling correlated stacks of image objects. 
The Stack class can be used alone, but is designed to work best from within other Dorado 
components such as Ceres.
'''


class Stack:
    '''
    The Stack class handles correlated stacks of image objects and their associated data.
    Stack exceeds simply storing stacks of images in an array or list by expanding upon the array
    object to hold additional data about the image stack such as a WCS object, the enclosed target,
    and the correlated calibration frames.

    Attributes
    ----------

    data: CCDdata array
        Array of CCDdata images.

    flat: CCDdata
        Flatfield calibration frame for data stack.

    filter: str
        Name of filter data was collected in.

    times: array or list-like
        Time for each image as an 'astropy.Time' object.

    calibrated: Boolean
        Whether the data is calibrated or not. Default is 'None'
        
    aligned: Boolean
        Whether the data is aligned or not. Default is 'None'
    target: Zellars object

    alignTo: int


    '''
    def __init__(self, data, flat = None, filter = '', times = [], calibrated = None, aligned = None, target = None, alignTo = 0):
        self.data = data
        self.flat = flat
        self.filter = filter
        self.length = len(data)
        self.calibrated = calibrated
        self.aligned = aligned

        self.target = target
        self.target_info = {} # dictionary of values, call simbad?

        self.times = times 
        self.wcs = None
        self.alignTo = alignTo
        self.solved = None
        # include things like flux uncertainty etc.
        # save data
        # if flat is none, call clippy for help

        if self.filter == '':
            try:
                self.filter = data[0].header['filter']
            except:
                self.filter = ''
        
        if self.times == []:
            try: 
                self.get_times()
                print('Times set')
            except:
                self.times = []
     
    def get_times(self):
        times = []
        for im in self.data:
            times.append(Time(im.header['DATE-OBS'], format='fits'))
        self.times = times

    def get_target_info(self, target = None):
        if target != None:
            self.target = target
   