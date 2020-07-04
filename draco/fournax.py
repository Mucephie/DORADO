import numpy as np

def fournax(x, y, terms):
    # Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'.
    
#     # Fourier series domain
#     tau = (np.max(x) - np.min(x))
    # Compute real valued Fourier transform
    f = np.fft.rfft(y)
    # Null or zero coefficients above ammount of series "terms"
    # This corresponds to undesired high-frequency terms
    f[terms+1:] = 0
    # Collapse back into function space, result is smoothed Fourier curve
    Fcurve = np.fft.irfft(f)
    
    return Fcurve

def pour(y):
    # Fournax is an abbreviation of Fourier numerical astronomy extension, its name is a backronym styled to match the constellation 'fornax'.

    # Compute real valued Fourier transform
    f = np.fft.fft(y)
    # p = np.square(np.abs(f))
    p = np.square(np.abs(f))
    
    return p
    
def var_ephem(OBS, epoch, period):
    cycle_ref = (OBS-epoch)/period
    cycle = np.round(cycle_ref)
    OmC = OBS-(epoch+cycle*period)
    return cycle, OmC
# def powerranger
# def maxlight
# def properties

