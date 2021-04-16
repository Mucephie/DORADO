import warnings
warnings.filterwarnings('ignore')
from astropy.time import Time

class Stack:
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
   