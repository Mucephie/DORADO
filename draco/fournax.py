import numpy as np
'''
Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'.
'''

__all__ = ['four_curve', 'pour', 'var_ephem']

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
# def powerranger
# def maxlight
# def properties

