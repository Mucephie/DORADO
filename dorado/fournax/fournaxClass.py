import numpy as np
from .zellars import Zellars
'''
Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'. 

'''

__all__ = ['Fournax']



class Fournax(Zellars):
    '''
    The Fournax class extends the Zellars target class to provide a consistent simple, yet robust interface to targets with regular or semi-regular photometric variability for the purposes of lightcurve/timeseries fourier analysis.  

    Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'. 

    Attributes
    ----------

    name: str
        name of target in string format.

    epoch: float or astropy.Time
        Epoch for the targets ephemeris of which the theoretical times of extrema will be extrapolated from.

    period: float
        The period between extrema of interest corresponding to the ephemeris epoch given. 

    '''
    def __init__(self, name, epoch = None, period = None):  
        Zellars.__init__(name)

        self.toml =[]
        self.

        if (epoch == None) and (period == None):
            print('No ephemeris data given for target')

        elif (epoch != None) and (period != None):
            self.epoch = epoch
            self.period = period # Must have units (TODO:: consider falling back to default unit if none)
            ## should these values be combined into an ephemeris dictionary to accomodate for an observational period to be determined later?

        elif (epoch != None):
            self.epoch = epoch 
            print('Period data not given.')
        elif (period != None):
            self.period = period
            print('Epoch data not given.')

        else:
            print('Unknown ephemeris error encountered. Please report via Dorado Github issues page.')
            print('Period given: ', period)
            print('Epoch given: ', epoch)

        
    def var_ephem(self, OBS):
        '''
        var_ephem takes observed time(s) of max light for a repeating variable star and ephemeris data and returns O-C values as well as the corresponding cycle.
        Parameters
        ----------
        OBS: float, list, or array
                Observed time(s) of max light
        epoch: float
                Epoch of ephemeris
        period: float
                Period of ephemeris 
        Returns
        -------
        cycle: float, list, or array
                The cycle corresponding to the time(s) of max light
        OmC: float, list, or array
                The O-C value for the time(s) of max light
        '''
        
        cycle_ref = (OBS - self.epoch) / self.period
        cycle = np.round(cycle_ref)
        OmC = OBS - (self.epoch + cycle * self.period)

        return cycle, OmC
        
    
    def superfit(self, filter, terms, s):
        '''
        superfit takes a raw timeseries and performs a multistep smoothed fit of the data
        which includes a spline fit tailored with the curvature of knots from a spline
        fit of a Fourier fit. The result is then expanded and then convoluted with a
        'blackman' window function.
        Parameters
        ----------
        x, y: array
                timeseries arrays to be fit
        terms: int
                terms to retain in the fourier fit.
        s: int
                new data array size 
        Returns
        -------
        X, Y: array
                The superfit timeseries arrays.
        '''

        ## TODO:: accomodate flux uncertainties
        ## TODO:: Zero mean of signal.
        ## TODO:: check if flux is photons or photons per second

        y = self.ts[self.filters[filter]]['flux']
        x = self.ts[self.filters[filter]]['time'] ## TODO:: convert to float friendly format like mjd

        ########################################################################
        # # Fourier fit the data to model curvature
        # f = np.fft.rfft(y)
        # # Null or zero coefficients above ammount of series "terms"
        # # This corresponds to undesired high-frequency terms
        # f[terms+1:] = 0
        # # Collapse back into function space, result is smoothed Fourier curve
        # F = np.fft.irfft(f)
        ########################################################################

        F = np.fft.irfft((np.fft.rfft(y)[terms+1:] = 0))
        # Create a spline fit of the fourier fit to extract knots
        tispl = InterpolatedUnivariateSpline(x[:len(F)], F, k = 5)
        # feed knots into spline of raw data
        LSQspl = LSQUnivariateSpline(x, y, tispl.get_knots()[1:-1]) 


        X = np.linspace(np.min(x), np.max(x), s)
        Y = self.smooth(LSQspl(X), window = 'blackman')[5:-5]

        ## TODO:: internalize the results
        # return X, Y


    def smooth(self, x, window_len=11, window='hanning'):
        """smooth the data using a window with requested size.
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                flat window will produce a moving average smoothing.
        output:
        the smoothed signal
        example:
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        see also: 
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
        2017-07-13 (last modified), 2006-10-31 (created)
        Section author: Unknown[1], GaelVaroquaux, Unknown[142], Unknown[143], Unknown[144], Unknown[145], Unknown[146], Unknown[147], WesTurner, Christian Gagnon, clecocel
        """
        ## Should this be removed from the class and instead be relegated to being a utility? Can this utiity be located in Dorado dependencies?
        if x.ndim != 1:
                raise ValueError("smooth only accepts 1 dimension arrays.")

        if x.size < window_len:
                raise ValueError("Input vector needs to be bigger than window size.")


        if window_len<3:
                return x


        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError("Window is none of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


        s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]
        # print(len(s))
        if window == 'flat': # moving average
                w = np.ones(window_len,'d')
        else:
                w = eval('np.' + window +'(window_len)')

        y = np.convolve(w/w.sum(), s, mode='valid')

        return y