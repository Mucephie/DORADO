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
        # date, location, weather, timezone, camera
        

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
        series = self.data[filter]
        for i in range(len(series)):
            if (operator == '+'):
                series[i].data = series[i].data  + operand
            elif (operator == '-'):
                series[i].data = series[i].data - operand
            elif (operator == '/'):
                series[i].data = series[i].data  / operand
            elif (operator == '*'):
                series[i].data = series[i].data  * operand
        
        self.data[filter] = series

    # save to wrk



# needs a way to save the stacks as one object and call them by name








class Stack:
    def __init__(self, data, filter = '', calibrated = None, aligned = None):
        self.data = data
        self.filter = filter
        self.length = len(data)
        self.calibrated = calibrated
        self.aligned = aligned
        # include things like flux uncertainty etc.
        # include times
        # include wcs
        # include target
        # include CCD temp

        if filter == '':
            try:
                filter = data[0].header['filter']
            except:
                filter = ''
        