import numpy as np
import ccdproc
'''
Ceres is the handler of series in Dorado,
'''

# __all__ = [Ceres]

class Ceres:

    def __init__(self):
        # calibrate
        # metadata
        self.filters = {}
        self.data = []
        self.bias = None
        self.flats = {}
        self.darks = {}
        self.time = {}
        # date, location, weather, timezone, camera, observer array
        # refering instance of clippy to call and save later?
        # or call clippy directly and feed it a ceres object
        

    def add_stack(self, stack):
        # eventually stacks themelves should have some metadata 
        # to denote stuff like calibration status
        self.filters[stack.filter] = len(self.data)
        self.data.append(stack)
        # this should also extract the time strings

    def rem_stack(self, filter):
        del self.data[self.filters[filter]]
        # delete time strings

    def calibrate(self, filter):
        # for bla in series: add bias corrected = True to header
        flat = self.flats[filter]
        bias = self.bias
        proc = self.data[filter]
        for p in range(len(proc)):
            proc[p] = ccdproc.ccd_process(proc[p], master_bias = bias, master_flat = flat)
        self.data[filter] = proc
    
    def imarith(self, filter, operator, operand):
        # mod to check datatype using type()
        # mod to remove im_count and make possible to use single image
        # mod to accomodate CCDdata object
        series = self.data[self.filters[filter]]
        for i in range(len(series)):
            if (operator == '+'):
                series[i].data = series[i].data  + operand
            elif (operator == '-'):
                series[i].data = series[i].data - operand
            elif (operator == '/'):
                series[i].data = series[i].data  / operand
            elif (operator == '*'):
                series[i].data = series[i].data  * operand
        
        self.data[self.filters[filter]] = series

    # save to wrk



class Stack:
    def __init__(self, data, filter = '', times = [], calibrated = None, aligned = None, target = ''):
        self.data = data
        self.filter = filter
        self.length = len(data)
        self.calibrated = calibrated
        self.aligned = aligned
        self.target = target
        self.target_info = {} # dictionary of values
        self.times = times # put together a check for if filled, if not, try to find it.

        # include things like flux uncertainty etc.
        # include wcs


        if filter == '':
            try:
                filter = data[0].header['filter']
            except:
                filter = ''
        