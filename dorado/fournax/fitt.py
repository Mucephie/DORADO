import numpy as np
from scipy.interpolate import LSQUnivariateSpline, InterpolatedUnivariateSpline, UnivariateSpline
'''
Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'. 

A major future update will bring a more powerful curve fitting function to this sub-package.
'''

# __all__ = ['four_curve', 'pour', 'var_ephem']
# class Fournax():

def four_curve(y, terms):
        """
        Four_curve takes a timeseries and 'smoothes' out higher order terms in frequency space and returns the smoothed values.

        Parameters
        ----------
        y: list or array
                y values of a timeseries. Time is x.
        terms: int
                terms to retain in fourier series.
        Returns
        -------
        Fcurve: list or array
                Fourier smoothed y values
        """
        # Fourier series domain
        # Compute real valued Fourier transform
        f = np.fft.rfft(y)
        # Null or zero coefficients above ammount of series "terms"
        # This corresponds to undesired high-frequency terms
        f[terms+1:] = 0
        # Collapse back into function space, result is smoothed Fourier curve
        Fcurve = np.fft.irfft(f)
        
        return Fcurve

def pour(y):
        '''
        pour takes a timeseries and computes the power in frequency space.

        Parameters
        ----------
        y: list or array
                y values of a timeseries. Time is x.

        Returns
        -------
        p: list or array
                power of y

        '''
        # Compute real valued Fourier transform
        f = np.fft.fft(y)
        # p = np.square(np.abs(f))
        p = np.square(np.abs(f))
        
        return p

def var_ephem(OBS, epoch, period):
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
        cycle_ref = (OBS-epoch)/period
        cycle = np.round(cycle_ref)
        OmC = OBS-(epoch+cycle*period)
        return cycle, OmC

def superfit(x, y, terms, s):
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

        # Fourier fit the data to model curvature
        Fcurve = four_curve(y, terms)
        # Create a spline fit of the fourier fit to extract knots
        tispl = InterpolatedUnivariateSpline(x[:len(Fcurve)], Fcurve, k = 5)
        # feed knots into spline of raw data
        LSQspl = LSQUnivariateSpline(x, y, tispl.get_knots()[1:-1]) 
        X = np.linspace(np.min(x), np.max(x), s)

        Y = smooth(LSQspl(X), window = 'blackman')

        Y = Y[5:-5]

        return X, Y

def smooth(x,window_len=11,window='hanning'):
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

# def maxlight
# def properties

