import numpy as np
import ccdproc
'''
Ceres is the handler of series in Dorado,
'''

__all__ = []

class ceres:

    def __init__(self):
        # calibrate
        # metadata
        self.data = {}
        self.bias = False
        self.flats = {}
        self.darks = {}
        self.time = {}
        

    def add_stack(self, stack, filter):
        # eventually stacks themelves should have some metadata 
        # to denote stuff like calibration status
        self.data[filter] = stack
        # this should also extract the time strings

    def rem_stack(self, filter):
        del self.data[filter]
        # delete time strings

    def calibrate(self, filter):
        # for bla in series: add bias corrected = True to header
        flat = self.flats[filter]
        bias = self.bias
        proc = self.data[filter]
        for p in range(len(proc)):
            proc[p] = ccdproc.ccd_process(proc[p], master_bias = bias, master_flat = flat)
        self.data[filter] = proc












# class stack:
#     def __init__(self, data, filter):
#         self.length = len(data)
#         # include things like flux uncertainty etc.